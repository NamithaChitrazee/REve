import subprocess

# launch pickEvent
result = subprocess.run(['pickEvent', '-e','-v','mcs.sophie.MDS2aTriggered.v0.art',str(args.run),'/',str(args.subrun),'/',str(args.run)], capture_output=True, text=True)
print(result.stdout)

