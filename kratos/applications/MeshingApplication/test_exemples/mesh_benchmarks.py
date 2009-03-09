import os

def Run():
	Msg = ""
	Text = "== MESHING_APPL ==========\n"

	# dam2d

	Text += "remesh: "
	os.chdir("adaptive_mesher2d.gid")	
	
	import remesh_benchmark.py
	Msg = remesh_benchmark.Run()
	
	if (Msg == True):
		Text += "OK\n"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"

	os.chdir("..")

	# Add other examples here

	return Text
