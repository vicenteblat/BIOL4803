import argparse
parser = argparse.ArgumentParser(description = 'This is a sample program')
parser.add_argument("loops", type=int, help="Enter the number of times you want a loop run")
parser.add_argument("message", help = "Enter a message")
args = parser.parse_args()
for i in range(0,args.loops):
    print("This is the " + str(i+1) + " time we ran the loop " + args.message)
