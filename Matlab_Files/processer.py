import matlab

def main():
    eng = matlab.engine.start_matlab("-nodisplay")
    eng.processor(nargout=0)

if __name__=='__main__':
    main()