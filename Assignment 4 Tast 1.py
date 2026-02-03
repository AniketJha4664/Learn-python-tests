try:
    a=open("sample.txt", 'rt')
    print("Reading file contents.")
    b=a.readline()
    c=a.readline()
    print(f"Line 1: {b}")
    print(f"Line 2: {c}")
except FileNotFoundError:
    print("Error: The file 'sample.txt' was not found.")


