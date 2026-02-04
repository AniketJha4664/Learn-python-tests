
score_dict={"Ritam":67, 'Aman':69, 'Vinayak':77, 'Abhijeet':62, 'Alice':85, 'Lily':12, 'Rashika':56, 'Riya':99, 'Sakshi':100}

input_name=input("Enter student's name: ")

presence_variable=input_name in score_dict.keys()
if presence_variable==True:
    result=score_dict.get(input_name)
    print(f"Students makrs is {result}")
else:
    print("Student not found.")

