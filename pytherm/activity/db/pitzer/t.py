
df = open(r'C:\Users\s215013\Documents\vscode\термодинамика\pytherm\pytherm\activity\db\pitzer\b.txt', 'r')
bufer = []
for line in df:
    b = line
    while '  ' in b:
        b = b.replace('  ', ' ')
    if b[0] == ' ':
        b = b[1:]
    if b[0] == '#':
        continue
    if b[-1:] == '\n':
        b = b[:-1]
    if b[-1] == ' ':
        b = b[:-1]
    bufer.append(b.split(' '))
df.close()

for line in bufer:
    if '+' in line[2]:
        line[2], line[1] = line[1], line[2]
    if '+' in line[1]:
        line[1], line[0] = line[0], line[1]
    if '+' in line[2]:
        line[2], line[1] = line[1], line[2]

bufer = sorted(bufer, key=lambda x: x[0])
bufer2 = []
for line in bufer:
    if '#' in line:
        i = line.index('#')
        b = ' '.join(line[i:])
        line = line[:i]
        line.append(b)
    bufer2.append(line)
    print(line)
bufer = bufer2

of = open(r'C:\Users\s215013\Documents\vscode\термодинамика\pytherm\pytherm\activity\db\pitzer\out.txt', 'w')
for line in bufer:
    b = ' '.join(line[3:])
    out = '\t'
    for i in line:
        out += f'{i:<13}'
    # of.write(f'\t{line[0]:<15}{line[1]:<15}{line[2]:<15}{b}\n')
    out += '\n'
    of.write(out)