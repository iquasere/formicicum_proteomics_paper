from subprocess import run, Popen, PIPE


def run_command(bash_command, print_message=True):
    if print_message:
        print(bash_command)
    run(bash_command.split(), check=True)


def run_pipe_command(bash_command, output='', mode='w', print_message=True):
    if print_message:
        print(bash_command)
    if output == '':
        Popen(bash_command, stdin=PIPE, shell=True).communicate()
    elif output == 'PIPE':
        return Popen(bash_command, stdin=PIPE, shell=True, stdout=PIPE).communicate()[0].decode('utf8')
    else:
        with open(output, mode) as output_file:
            Popen(bash_command, stdin=PIPE, shell=True, stdout=output_file).communicate()


def parse_fasta(file):
    print(f'Parsing {file}')
    with open(file) as f:
        lines = [line.rstrip('\n') for line in f]
    i = 0
    sequences = {}
    while i < len(lines):
        if lines[i].startswith('>'):
            name = lines[i][1:]
            sequences[name] = ''
            i += 1
            while i < len(lines) and not lines[i].startswith('>'):
                sequences[name] += lines[i]
                i += 1
    return sequences
