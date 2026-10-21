import json
import socket
import struct

import numpy as np

import KratosMultiphysics as KMP

def Factory(settings, Model):
    if not isinstance(settings, KMP.Parameters):
        raise TypeError("Expected a Kratos Parameters object")
    return SocketCommunicationProcess(Model, settings["Parameters"])

class SocketCommunicationProcess(KMP.Process):
    """ This process allows to send and receive data from a socket connection. 
        It can be used to exchange data with external programs or scripts during the simulation. 
    """
    
    # Defines the list of hooks where the recv and send actions can be executed.
    # ExecuteInitialize and ExecuteFinalize are intentionally left ousdie the list of hooks
    # because they are reserved for the handshake and cleanup, respectively.
    _hooks = (
        "Execute", 
        "ExecuteBeforeSolutionLoop",
        "ExecuteInitializeSolutionStep", 
        "ExecuteFinalizeSolutionStep",
        "ExecuteBeforeOutputStep", 
        "ExecuteAfterOutputStep", 
    )

    @classmethod
    def GetDefaultParameters(cls):
        return """{
            "model_part_name": "",
            "host": "127.0.0.1",
            "port": 0,
            "timeout": 30.0,
            "recv": {
                "variable_name": "",
                "execute_at": "ExecuteInitializeSolutionStep",
                "historical": true
            },
            "send": {
                "variable_name": "",
                "execute_at": "ExecuteFinalizeSolutionStep",
                "historical": true
            }
        }"""

    # For the future
    @classmethod
    def GetSchema(cls):
        return """{
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "type": "object",
            "properties": {
                "model_part_name": {
                    "type": "string",
                    "default": ""
                },
                "host": {
                    "type": "string",
                    "default": "127.0.0.1"
                },
                "port": {
                    "type": "integer",
                    "minimum": 0,
                    "maximum": 65535,
                    "default": 0
                },
                "timeout": {
                    "type": "number",
                    "minimum": 0,
                    "default": 30.0
                },
                "recv": {
                    "$ref": "#/$defs/communication",
                    "properties": {
                        "execute_at": {
                        "default": "ExecuteInitializeSolutionStep"
                        }
                    }
                    },
                "send": {
                    "$ref": "#/$defs/communication",
                    "properties": {
                        "execute_at": {
                        "default": "ExecuteFinalizeSolutionStep"
                        }
                    }
                }
            },
            "$defs": {
                "communication": {
                    "type": "object",
                    "properties": {
                        "variable_name": {
                            "type": "string",
                            "default": ""
                        },
                        "execute_at": {
                            "type": "string"
                        },
                        "historical": {
                            "type": "boolean",
                            "default": true
                        }
                    },
                    "required": [
                        "variable_name",
                        "execute_at",
                        "historical"
                    ],
                    "additionalProperties": false
                }
            },
            "required": [
                "model_part_name",
                "host",
                "port",
                "timeout"
            ],
            "additionalProperties": false
        }"""

    def __init__(self, model, settings):
        super().__init__()

        defaults = KMP.Parameters(SocketCommunicationProcess.GetDefaultParameters())

        settings.ValidateAndAssignDefaults(defaults)
        
        for direction in ("recv", "send"):
            settings[direction].ValidateAndAssignDefaults(defaults[direction])
        
        self.model_part = model[settings["model_part_name"].GetString()]
        self.host = settings["host"].GetString()
        self.port = settings["port"].GetInt()
        self.timeout = settings["timeout"].GetDouble()

        if not self.host or not 0 < self.port < 65536 or self.timeout <= 0:
            raise ValueError("A host, TCP port (1-65535), and positive timeout are required")

        self.actions = {}

        for direction in ("recv", "send"):
            action = settings[direction]
            name = action["variable_name"].GetString()    
            hook = action["execute_at"].GetString()
            historical = action["historical"].GetBool()

            if not name:
                raise ValueError(f"A variable name or COORDINATES must be specified.")

            if hook not in self._hooks:
                raise ValueError(f"Unsupported execute_at '{hook}' for '{direction}'")
            
            if name == "COORDINATES":
                variable = None
                components = 3
            else:
                variable_type = KMP.KratosGlobals.GetVariableType(name)
                
                if variable_type not in ("Double", "Array"):
                    raise ValueError(f"Only Double and Array<3, double> variables are supported: {name}")
                
                variable = KMP.KratosGlobals.GetVariable(name)
                components = 1 if variable_type == "Double" else 3
            
            self.actions[direction] = (name, variable, components, hook, historical)
        
        if not self.actions:
            raise ValueError("Configure at least one of receive or send")

        self.connection = None
        self.node_ids = None
        self.adaptors = {}

    ################################################################
    #### Reserved hooks                                         ####
    ################################################################

    def ExecuteInitialize(self):
        self.node_ids = [node.Id for node in self.model_part.Nodes]
        
        # Initialize the adaptors for the type of data we want to communicate
        for direction, (name, variable, _, _, historical) in self.actions.items():
            if name == "COORDINATES":
                adaptor = KMP.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, KMP.Configuration.Current)
            elif historical:
                adaptor = KMP.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, variable, 0)
            else:
                adaptor = KMP.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, variable)
            self.adaptors[direction] = adaptor

        # Create the connection with the device
        self.connection = socket.create_connection((self.host, self.port), self.timeout)
        
        try:
            self._handshake()
            
        except Exception:
            self.connection.close()
            self.connection = None
            raise

    def ExecuteFinalize(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None

    ################################################################
    #### Available hooks                                        ####
    ################################################################

    def Execute(self):
        self._execute_actions("Execute")

    def ExecuteBeforeSolutionLoop(self):
        self._execute_actions("ExecuteBeforeSolutionLoop")

    def ExecuteInitializeSolutionStep(self):
        self._execute_actions("ExecuteInitializeSolutionStep")

    def ExecuteFinalizeSolutionStep(self):
        self._execute_actions("ExecuteFinalizeSolutionStep")

    def ExecuteBeforeOutputStep(self):
        self._execute_actions("ExecuteBeforeOutputStep")

    def ExecuteAfterOutputStep(self):
        self._execute_actions("ExecuteAfterOutputStep")

    ################################################################
    #### Impl                                                   ####
    ################################################################

    def _execute_actions(self, hook):
        for direction in ("recv", "send"):
            if direction not in self.actions or self.actions[direction][3] != hook:
                continue
            
            if self.connection is None:
                raise RuntimeError("SocketDataProcess must be initialized before data exchange")
            
            if [node.Id for node in self.model_part.Nodes] != self.node_ids:
                raise RuntimeError("Model part node IDs changed after the socket handshake")
            
            name, _, components, _, _ = self.actions[direction]
            
            shape = (len(self.node_ids), components)
            
            header = {
                "type": "data", 
                "variable_name": name,
                "shape": list(shape), 
                "dtype": "<f8"
            }
            
            adaptor = self.adaptors[direction]
            
            if direction == "recv":
                if self._receive_header() != header:
                    raise ValueError("Unexpected socket data header for " + name)

                data = np.frombuffer(self._receive_exact(8 * shape[0] * shape[1]), dtype="<f8").reshape(shape)
                adaptor.data = data[:, 0] if components == 1 else data
                adaptor.StoreData()

            if direction == "send":
                adaptor.CollectData()
                data = np.asarray(adaptor.data, dtype="<f8").reshape(shape)
                self._send_header(header)
                self.connection.sendall(data.tobytes(order="C"))

    ################################################################
    #### Comm Layer                                             ####
    ################################################################

    def _handshake(self):
        """ Performs the initial handshake for the communication and exchanges 
            simulation data in case there is any """

        header = {
            "type": "handshake", 
            "version": 1,
            "number_of_nodes": len(self.node_ids)
        }

        self._send_header(header)
        reply = self._receive_header()
        
        # TODO: Compare something more robust, like a checksum.
        if reply != header:
            raise ValueError("Socket peer handshake does not match the model part nodes")

        # TODO: Exchange initial data...

    def _send_header(self, header):
        encoded = json.dumps(header, separators=(",", ":")).encode("utf-8")
        self.connection.sendall(struct.pack("!I", len(encoded)) + encoded)

    def _receive_header(self):
        size = struct.unpack("!I", self._receive_exact(4))[0]    
        return json.loads(self._receive_exact(size).decode("utf-8"))

    def _receive_exact(self, size):
        chunks = bytearray(size)
        view = memoryview(chunks)
        received = 0
        
        # Recieve the info by chunks
        while received < size:
            count = self.connection.recv_into(view[received:])
            
            if count == 0:
                raise ConnectionError("Socket closed during data exchange")
            received += count

        return chunks
