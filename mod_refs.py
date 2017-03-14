import sys

if __name__ == "__main__":

    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
            if line.startswith(">"):
                tokens = line.strip().strip(">").split("|")
                tmp = tokens[-1]
                tokens[-1] = tokens[0]
                tokens[0] = tmp
                print ">" + "|".join(tokens)
            else:
                print line.strip()
