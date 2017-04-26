class Utility(object):
    @staticmethod
    def writefile(pathandname, sequence):
        f = open(name=pathandname, mode='w')
        f.write(sequence)
        f.close()

    @staticmethod
    def readfile(pathandname):
        sequence = []
        with open(name=pathandname, mode='r') as f:
            for line in f:
                sequence.append(line)

        f.close()
        return ''.join(sequence)
