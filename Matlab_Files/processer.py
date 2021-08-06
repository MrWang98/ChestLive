import os

import matlab.engine

def main():
    eng = matlab.engine.start_matlab("-nodisplay")
    eng.processor(nargout=0)

if __name__=='__main__':
    print(os.path.abspath(os.curdir))
    main()
