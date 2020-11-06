#--------------------------------------------------------------------------
#	ASSIGNMENT 1 - Introduction to Python
#	DUE: 08/29/2019 at 9:00pm
#--------------------------------------------------------------------------

# 	1A.	In the string class, find two string methods, METHOD1 and METHOD2
#		that will return True if the string can be converted into an
#		integer. (5 pts)


testCase1 = '23432'
testCase1.isnumeric() # should return TRUE
testCase1.isdigit() # should return TRUE


testCase2 = '23423.234'
testCase2.isnumeric() # should return FALSE
testCase2.isdigit() # should return FALSE


testCase3 = 'PatrickMcGrath'
testCase3.isnumeric() # should return FALSE
testCase3.isdigit() # should return FALSE



# 	1B. Using Google, determine what the difference is between these two
#		methods. (5 points)
### The difference between these two methods is that isdigit() returns True
### if and only if all the characters in the string are digit characters, that
### is if they are a character from 0 to 9 inclusive. On the other hand,
### isnumeric() returns True if all the characters are numeric, which includes
### fractions for example. Eg. using the String "string = '\u00BD'" (which is
### the fraction 1/2) isnumeric() would return True but isnumeric() would
### return False.

#--------------------------------------------------------------------------

# 	2. 	Sometimes you want to pad a string so that it’s a certain length.
#		For example, let’s say you have a DNA sequence that you wanted to
#		be 9 bp long, adding ‘N’s for any bases that are missing. Find two
#		methods that would allow you to do the following: (5 pts each)


b = 'CGTTT'
b.ljust(9,'-') # should return 'CGTTT----'
b.rjust(9,'-') # should return '----CGTTT'


#--------------------------------------------------------------------------

#	3. 	Assume you have a float, 3454.4521. You want to split this number
#		into the integer and decimal portion of the number. (5 pts each)
#		A. 	In the first case, just use methods or magic methods of the
#			float class to accomplish this.
#		B. 	Convert the float to a string and use string methods to accomplish
#			this (the final numbers should be integers)


a = 3454.4521
a.__int__() # should return 3454
a.__divmod__(1)[1] # should return .4520999999999731


b = str(a)
int(b.split('.')[0]) # should return 3454
int(b.split('.')[1]) # should return 4521


#--------------------------------------------------------------------------

# 	4. 	Create a function that takes in three numbers and returns the sum
#		(10 pts)


a = 8
b = 29
c = 2019

def calculator(a,b,c):
    return(a + b + c)


#--------------------------------------------------------------------------

#	5. 	Create a function that takes in a string containing 10 words
#		separated by commas and returns a list of the words in alphabetical
#		order (10 pts)


a = 'python,string,this,is,a,test,of,the,georgia,tech'
def spliter(str):
    return(str.split(','))


#--------------------------------------------------------------------------

#	6. 	Let a = True
#		A.	What kind of object is a? (4 pts)
#		B. 	What kind of object is returned by a + a? (3 pts)
#		C. 	What type of object is returned by a * a? (3 pts)


### A. Boolean
### B. Integer
### C. Integer


#--------------------------------------------------------------------------

#	7. 	Assume you have a string of valid DNA sequence (all upper case).
#		DNA methylation enzymes can add methyl groups to DNA nucleotides
#		at defined sequences. However, there can be redundancy in the motif
#		recognized by the protein. For example, the CcrM methyltransferase
#		recognizes GANTC sites (where N can be any of the four nucleotides).
#		Create a single line that returns the number of CcrM methylation
#		sites within the string of DNA sequence (call the variable
#		DNA_sequnce). (10 pts)


DNA_sequence = 'AGAAGGCCCTAGAGTCCAAGACTCAGATTCAGAATCGACTC'
DNA_sequence.count('GAATC') + DNA_sequence.count('GAGTC') + DNA_sequence.count('GACTC') + DNA_sequence.count('GATTC')


#--------------------------------------------------------------------------

# 	8. 	You are asked to take a string of DNA sequence and return a string
#		of RNA sequence by replacing all the Ts with Us. Do this in two ways:
#		A. 	Use a single string method to accomplish this (5 pts)
#		B. 	Use a combination of string’s split and join methods to accomplish
#			this (5 pts)


DNA_sequence = 'AGAAGGCCCTAGAGTCCAAGACTCAGATTCAGAATCGACTC'


# 8A:
DNA_sequence.replace('T','U') # should return 'AGAAGGCCCUAGAGUCCAAGACUCAGAUUCAGAAUCGACUC'



# 8B:
mylist = DNA_sequence.split('T')
RNA_sequence = 'U'
RNA_sequence.join(mylist)


#--------------------------------------------------------------------------

#	9. 	These are 2.5 pts each:
#		A. 	Create a list of numbers between 1 and 15.
#		B. 	Use a list comprehension to convert the numbers to strings and
#			add a letter x in front of each
#		C. 	Use a string method to concatenate the list into a single
#			string separated by two periods
#		D. 	Chain all of these commands together into a single line that
#			accomplishes all three tasks


# 9A-C:
a = list(range(1,16))
b = ['x' + i.__str__() for i in a]
c = '..'.join(b)


# 9D:
'..'.join(['x' + i.__str__() for i in list(range(1,16))])


#--------------------------------------------------------------------------

# 	10. You are given a list of tasks by your professor (2.5 pts ea):
#		A.	Sort the list in alphabetical order
#		B. 	You are given a string representing an additional single task.
#			Add this task to the beginning of the list.
#		C. 	Add this task to the end of the list
#		D.	You are now given a new list of tasks. Add these tasks to the
#			list of tasks.


# 10A:
tasks = ['Read Lecture', 'Go to class', 'Work on homework']
tasks.sort()


# 10B:
new_task = 'Start computer'
tasks.insert(0,new_task)
tasks # should return ['Start computer', 'Go to class', 'Read Lecture', 'Work on homework']


# 10C:
tasks.pop(0)
tasks.append(new_task)
tasks # should return ['Go to class', 'Read Lecture', 'Work on homework','Start Computer']


# 10D:
new_tasks = ['Study', 'Check grade']
tasks.extend(new_tasks)
tasks # should return ['Go to class', 'Read Lecture', 'Work on homework', 'Start computer', 'Study', 'Check grade']
