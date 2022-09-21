import subprocess as sp


def runIPANEMAP(configfile):
    cmd = ["ipanemap", "--config", configfile]
    
    #[print(arg,end=" ") for arg in cmd]
    sp.run(cmd)
