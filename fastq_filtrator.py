import sys


def parse_arguments() -> dict:
    args = sys.argv
    arguments = {}
    options = ("--min_length", "--keep_filtered", "--gc_bounds", "--output_base_name")

    for idx, arg in enumerate(args):
        if arg == "--min_length":
            try:
                arguments[arg.lstrip("-")] = int(args[idx+1])
            except ValueError:
                print("Invalid value for --min_length. Must be integer.")
                sys.exit(1)
            break

    else:
        arguments["min_length"] = 0

    for idx, arg in enumerate(args):
        if arg == "--keep_filtered":
            arguments["keep_filtered"] = True
            break
    else:
        arguments["keep_filtered"] = False

    for idx, arg in enumerate(args):
        if arg == "--gc_bounds":
            arguments["gc_bounds"] = [0, 100]
            try:
                arguments["gc_bounds"][0] = int(args[idx+1])
                if args[idx+2].isdigit():
                    arguments["gc_bounds"][1] = int(args[idx+2])
            except ValueError:
                print("Invalid value for --gc_bounds. Must be integers from 0 to 100.")
                sys.exit(1)
            # Only if --gc_bounds is last argument in command and only 1 boundary specified.
            except IndexError:
                pass
            break
    else:
        arguments["gc_bounds"] = [0, 100]

    for idx, arg in enumerate(args):
        if arg.endswith(".fastq"):
            arguments["file_name"] = arg
            break
    else:
        print("Invalid file format. File must have .fastq extention")
        sys.exit(1)

    for idx, arg in enumerate(args):
        if arg == "--output_base_name":
            arguments["output_base_name"] = args[idx+1]
            break
    else:
        arguments["output_base_name"] = arguments["file_name"].split(".")[0]

    for option in sys.argv[1:]:
        if option.startswith("--") and option not in options and option != arguments["file_name"]:
            print(f"Unknown option {option}")
            sys.exit(1)

    return arguments


def fastq_filtrator(agruments: dict):
    out_passed = open(f'{arguments["output_base_name"]}_passed.fastq', 'w')
    if agruments['keep_filtered']:
        out_failed = open(f'{arguments["output_base_name"]}_failed.fastq', 'w')

    with open(agruments['file_name']) as fastq_file:
        while True:
            seq_name = fastq_file.readline()
            seq = fastq_file.readline()
            seq_length = len(seq[:-1])
            comment = fastq_file.readline()
            quality = fastq_file.readline()

            if not seq_name:
                break

            if seq_length < arguments['min_length']:
                if arguments['keep_filtered']:
                    out_failed.write(seq_name)
                    out_failed.write(seq)
                    out_failed.write(comment)
                    out_failed.write(quality)
                continue

            if arguments['gc_bounds'] != (0, 100):
                gc_quantity = seq.upper().count('C') + seq.upper().count('G')
                gc_content = gc_quantity / seq_length * 100
                if not (arguments['gc_bounds'][0] <= gc_content <= arguments['gc_bounds'][1]):
                    if arguments['keep_filtered']:
                        out_failed.write(seq_name)
                        out_failed.write(seq)
                        out_failed.write(comment)
                        out_failed.write(quality)
                    continue

            out_passed.write(seq_name)
            out_passed.write(seq)
            out_passed.write(comment)
            out_passed.write(quality)

    out_passed.close()
    if agruments['keep_filtered']:
        out_failed.close()

arguments = parse_arguments()

if not arguments['output_base_name']:
    print('''Вы указали аргумент --output_base_name, но не указали для него значение.\n
           Вы уверены, что хотите продолжить?''', 'y/n')
    f = input().lower()
    if f == 'y' or f == 'yes':
        fastq_filtrator(arguments)
    else:
        sys.exit(1)

if arguments['min_length'] == 0 and arguments['gc_bounds'] == (0, 100):
    print('''При данном наборе параметров все последовательности пройдут фильтрацию.\n
           Вы действительно хотите продолжить?''', 'y/n')
    f = input().lower()
    if f == 'y' or f == 'yes':
        fastq_filtrator(arguments)
    else:
        sys.exit(1)

if arguments['gc_bounds'][0] > arguments['gc_bounds'][1]:
    print('Нижняя граница --gc_bounds больше верхней, все последовательности не пройдут фильтрацию.')
    sys.exit(1)

fastq_filtrator(arguments)
