import sys

input = sys.argv[1]
output = sys.argv[2]

print(input)
print(output)
error_count = 0

with open(input,'r') as inp:
    out = open(output, 'w')
    print("open")
    for line in inp:
        linesplit = line.rstrip("\n").split('\t')
        if len(linesplit) < 3:
            continue
        try:
            if int(linesplit[1]) < int(linesplit[2]):
                out.write(line)
        except Exception as e:
            print(e)
            error_count += 1
            print(error_count)

    out.close()