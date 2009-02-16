import os

def Run():
	Msg = ""
	Text = "== PFEM ==========\n"

	# dam2d

	Text += "dam_2d: "
	os.chdir("dam2d.gid")	
	
	import dam2d_benchmark
	Msg = dam2d_benchmark.Run()
	
	if (Msg == True):
		Text += "OK\n"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"

	os.chdir("..")

	# Add other examples here

	return Text
