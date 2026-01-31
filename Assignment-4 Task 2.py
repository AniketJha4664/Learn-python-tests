with open("output.txt", 'wt') as a:
    input_1=input("Enter text to write to file: ")
    a.write(input_1)
    print("Data successfully written to output.txt")



with open("output.txt", 'at') as x:
    input_2=input("Enter additional text to append: ")
    x.write(f"\n{input_2}")
    print("Data successfully appended.")


with open("output.txt", 'rt') as r:
    line1=r.readline()
    print(line1)
    line2=r.readline()
    print(line2)

