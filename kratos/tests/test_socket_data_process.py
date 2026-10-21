import json
import socket
import struct

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics import socket_data_process

# This methods emulate a 3rd party code connecting towards the SocketDataProcess

def _read_exact(connection, size):
    data = bytearray(size)
    view = memoryview(data)
    received = 0
    while received < size:
        count = connection.recv_into(view[received:])
        if count == 0:
            raise ConnectionError("Peer closed the test connection")
        received += count
    return data


def _read_header(connection):
    size = struct.unpack("!I", _read_exact(connection, 4))[0]
    return json.loads(_read_exact(connection, size).decode("utf-8"))


def _send_header(connection, header):
    encoded = json.dumps(header).encode("utf-8")
    connection.sendall(struct.pack("!I", len(encoded)) + encoded)


def _peer(listener, received_name, received_values, sent_name, sent_shape):
    with listener.accept()[0] as connection:
        connection.settimeout(5.0)
        handshake = _read_header(connection)
        if handshake != {"type": "handshake", "version": 1, "number_of_nodes": 2}:
            raise AssertionError("Unexpected handshake: {}".format(handshake))
        _send_header(connection, handshake)

        values = np.asarray(received_values, dtype="<f8")
        _send_header(connection, {
            "type": "data", "variable_name": received_name,
            "shape": list(values.shape), "dtype": "<f8"
        })
        connection.sendall(values.tobytes(order="C"))

        header = _read_header(connection)
        expected_header = {
            "type": "data", "variable_name": sent_name,
            "shape": list(sent_shape), "dtype": "<f8"
        }
        if header != expected_header:
            raise AssertionError("Unexpected data header: {}".format(header))
        data = np.frombuffer(
            _read_exact(connection, 8 * np.prod(sent_shape)), dtype="<f8"
        ).reshape(sent_shape)
        return data.copy()


class TestSocketDataProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.model_part.CreateNewNode(2, 1.0, 2.0, 3.0)
        self.model_part.CreateNewNode(5, 4.0, 5.0, 6.0)

    # This methos sets up a SocketDataProcess and a "fake" extenral code trying to interact with it 
    def _exchange(self, send, recv, received_values, hooks):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as listener:
            
            listener.bind(("127.0.0.1", 0))
            listener.listen(1)
            listener.settimeout(5.0)

            settings = KM.Parameters(json.dumps({
                "Parameters": {
                    "model_part_name": "Main",
                    "host": "127.0.0.1",
                    "port": listener.getsockname()[1],
                    "timeout": 5.0,
                    "recv": recv,
                    "send": send
                }
            }))
            
            with ThreadPoolExecutor(max_workers=1) as executor:
                result = executor.submit(
                    _peer, 
                    listener, 
                    recv["variable_name"], 
                    received_values,
                    send["variable_name"],
                    (len(self.model_part.Nodes), 1 if send["variable_name"] == "PRESSURE" else 3)
                )

                process = socket_data_process.Factory(settings, self.model)
                
                try:
                    process.ExecuteInitialize()
                    for hook in hooks:
                        getattr(process, hook)()
                finally:
                    process.ExecuteFinalize()
                return result.result(timeout=5.0)

    def test_historical_scalar_receive_and_vector_send(self):
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.VELOCITY, 0, [node.Id, node.Id + 1, node.Id + 2])

        sent = self._exchange(
            # Send
            {
                "variable_name": "VELOCITY",
                 "execute_at": "ExecuteFinalizeSolutionStep",
                 "historical": True
            },
            # Recv
            {   
                "variable_name": "PRESSURE",
                 "execute_at": "ExecuteInitializeSolutionStep",
                 "historical": True
            },
            # Data
            [
                [7.5], [8.5]
            ],
            # Hooks
            (
                "ExecuteInitializeSolutionStep", 
                "ExecuteFinalizeSolutionStep"
            )
        )

        self.assertEqual(
            [node.GetSolutionStepValue(KM.PRESSURE) for node in self.model_part.Nodes],
            [7.5, 8.5]
        )
        np.testing.assert_allclose(sent, [[2.0, 3.0, 4.0], [5.0, 6.0, 7.0]])

    def test_coordinates_receive_and_nonhistorical_scalar_send(self):
        for node in self.model_part.Nodes: node.SetValue(KM.PRESSURE, node.Id * 2.0)

        sent = self._exchange(
            # Send
            {
                "variable_name": "PRESSURE",
                "execute_at": "ExecuteFinalizeSolutionStep", 
                "historical": False
            },
            # Recv
            {
                "variable_name": "COORDINATES", 
                "execute_at": "ExecuteInitializeSolutionStep",
                "historical": False
            },
            # Data
            [
                [10.0, 11.0, 12.0], 
                [13.0, 14.0, 15.0]
            ],
            # Hooks
            (
                "ExecuteInitializeSolutionStep", 
                "ExecuteFinalizeSolutionStep"
            )
        )

        np.testing.assert_allclose(
            [[node.X, node.Y, node.Z] for node in self.model_part.Nodes],
            [[10.0, 11.0, 12.0], [13.0, 14.0, 15.0]]
        )
        np.testing.assert_allclose(sent, [[4.0], [10.0]])


if __name__ == "__main__":
    KratosUnittest.main()
