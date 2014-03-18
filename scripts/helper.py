
def run_command(command):
    import subprocess

    runCmd = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output,errors = runCmd.communicate()
    print "\n" + command + "\n"
    print output
    print errors
