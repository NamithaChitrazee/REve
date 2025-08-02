import subprocess
import argparse

def main():
  # Step 1: launch pickEvent
  result = subprocess.run(['pickEvent', '-e','-v',str(args.dataset),str(args.run),'/',str(args.subrun),'/',str(args.run)], capture_output=True, text=True)
  print(result.stdout)
  
  # Step 2: TODO intermediate step to copy the selected file to a .txt
  
  # Step 3: launch the display
  launch_display = subprocess.run(['mu2e','-c','EventDisplay/examples/nominal_exampl.fcl','-S','file.txt' ], capture_output=True, text=True)


if __name__ == "__main__":
    # list of input arguments, defaults should be overridden
    parser = argparse.ArgumentParser(description='command arguments', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--dataset", type=str, required=True, help="mcs data set")
    parser.add_argument("--run", type=str, required=True, help="run number")
    parser.add_argument("--subrun", type=int, required=True,help="sub-run number")
    parser.add_argument("--event", type=int, required=True,help="event number")
    parser.add_argument("--verbose", default=1, help="verbose")
    args = parser.parse_args()
    (args) = parser.parse_args()

    # run main function
    main(args)
