# In order to make the MShePy module accessible either
#   - add the MIKE installation path/bin/x64 directory to PYTHONPATH
import MShePy as ms
setup_path = "../Data/3x3_Box/3x3_Box.she"

def main():
  ms.wm.runAll(setup_path)

if __name__ == "__main__":
  main()
