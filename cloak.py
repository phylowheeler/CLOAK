import os
import sys


'''
Name of File: cloak.py
Author: Chiragdeep Chatur (Masel Lab, University of Arizona)

Input Specifications:
    Create a directory containing multiple sequence alignments for the same alignment. Pass in the complete path to this directory as input for this program. 
Process:
    This program takes a set of multiple sequence alignments for the same set of sequences, and masks the positions for which there was not complete consensus.
    This program will also sort the species in each .fasta file as per the ID. Different alignment algorithms may product output in different orders, but this program
    addresses this issue by creating an intermediary directory called "ordered_fasta_files", which ensures that all .fasta files contain the species and their 
    corresponding sequences in the same order.  
Output Specifications:
    This program generates an output file named last_component+"_result.fasta" in the same directory from which it is run. This file contains the result, formatted in 
    fasta style. The last_component part indicates the input folder from which the result file was created. 
'''

alignmentFolder = sys.argv[1].rstrip("/")
#the alignmentFolder contains the full path to the directory in which multiple alignments of the same sequence are stored
#I am expecting this directory to contain fasta files

last_component = os.path.basename(alignmentFolder)
#this contains the last component of the input folder
#we use the last component to name the ordered_files folder and the result fasta file

 
'''
The fasta files provided may be out of order - the species may not be ordered in the same way in all fasta files. This function fixes this problem, and returns a 
directory of fasta files in which all fasta files have the same order of species
'''
def order_the_alignments(directory):
    new_dir = "ordered_fasta_files_"+last_component
    parent_dir = "./"
    path = os.path.join(parent_dir, new_dir)
    os.mkdir(path)
    for path, dirs, files in os.walk(directory):
        numFile = 0
        for f in files:
            myDict = {}
            fileName = os.path.join(path, f)
            with open(fileName, "r") as myFile:
                intermediarylines = []
                for line in myFile:
                    intermediarylines.append(line)
                lines = []
                current_line = ""

                for line in intermediarylines:
                    if line[0] == ">":
                        if current_line:
                            lines.append(current_line)
                        lines.append(line)
                        current_line = ""
                    else:
                        current_line += line

                if current_line:
                    lines.append(current_line)

                for i in range(len(lines)):
                    line = lines[i]
                    if len(line) == 0:
                        continue
                    if line[0] == ">":
                        if "|" in line:
                            second_half_string = line.split("|")[1]
                            species_id = second_half_string.split("-")[0]
                            myDict[species_id] = [line, lines[i+1]]
                        else:
                            #do what needs to be done for current format
                            species_id = line[1:]
                            myDict[species_id] = [line, lines[i+1]]
                fasta_filename = new_dir+"/fasta"+str(numFile)+".fasta"
                numFile += 1
                with open(fasta_filename, "a") as f:
                    for key in sorted(myDict.keys()):
                        f.writelines(myDict[key][0])
                        f.writelines(myDict[key][1])
    return [new_dir, myDict]


'''
This function takes in the full path to the directory containing the ordered fasta files as an input, and returns a 3d array. 
The 3d array is structured as follows:
    Each 2d array in the 3d array represents a different alignment
    Each 1d array in each 2d array represents a line in the alignment
'''
def convertToLetters(directory):
    return_array = []
    for path, dirs, files in os.walk(directory):
        for f in files:
            fileName = os.path.join(path, f)
            with open(fileName, "r") as myFile:
                intermediary_two_d_array = []
                for line in myFile:
                    line = line.strip()
                    line = list(line)
                    intermediary_two_d_array.append(line)
                two_d_array = []
                current_line = []
                for line in intermediary_two_d_array:
                    if line[0] == ">":
                        if len(current_line) != 0:
                            two_d_array.append(current_line)
                        current_line = []
                    else:
                        for char in line:
                            current_line.append(char)
                if len(current_line) != 0:
                    two_d_array.append(current_line) 
                return_array.append(two_d_array)
                two_d_array = []
    return return_array


'''
For the purpose of determining consensus, we need to convert the 3d array of letters and dashes to a 3d array of numbers and dashes. 
The numbering system I am using is borrowed from MergeAlign
'''
def convert_to_nums(array):
    nums = []
    for x in range(0, len(array)):
        two_d_array = array[x]
        make_two_d_array = []
        for y in range(0, len(two_d_array)):
            row = two_d_array[y]
            i = 1
            row_array = []
            for z in range(0, len(row)):
                letter = row[z]
                if letter == "-":
                    row_array.append("-")
                else:
                    row_array.append(i)
                    i += 1
            make_two_d_array.append(row_array)
            row_array = []
        nums.append(make_two_d_array)
        make_two_d_array = []    
    return nums


'''
The purpose of this method is to create the final output that our consensus algorithm is supposed to yield. 
'''
def create_output(all_als):
    output = []
    first_alignment = all_als[0]

    for x in range(len(first_alignment[0])):
        toSearch = create_toSearch(first_alignment, x) #create the initial toSearch term
        indices_of_all_alignments = [] #contains arrays of indices for every alignment

        for k in range(0, len(all_als)):
            indices = create_indices(all_als, toSearch, k)
            indices_of_all_alignments.append(indices)

        divvy_is_needed_for_current_column = divvy_is_required(indices_of_all_alignments)

        if divvy_is_needed_for_current_column == False:
            output.append(toSearch) #if none of the indices arrays are divvy-able, append toSearch to output
        if divvy_is_needed_for_current_column == True:
            add_to_output_or_discard(indices_of_all_alignments, output, toSearch)
        toSearch = []
    return output


'''
This function recursively adds or discards potential divvied results to output. 
'''
def add_to_output_or_discard(indices_of_all_alignments, output, toSearch, seen = []):
    number_of_nums = 0
    for num in toSearch:
        if num != "-":
            number_of_nums += 1
    if number_of_nums == 0:
        return
    array_containing_dicts_of_divvies = get_divvied_results_as_array_of_dicts(indices_of_all_alignments, toSearch)
    three_d_array_of_divvies = divvies_as_3darray(array_containing_dicts_of_divvies) 
    for i in range(1, len(three_d_array_of_divvies)):
        divvies_from_specific_alignment = three_d_array_of_divvies[i]
        for divvy in divvies_from_specific_alignment:
            can_this_divvy_be_found = can_this_divvy_be_found_without_subsequent_divvies(divvy, three_d_array_of_divvies, i)
            if toSearch == divvy:
                continue
            if can_this_divvy_be_found == True:
                output.append(divvy)
                continue
            if can_this_divvy_be_found == False:
                add_to_output_or_discard(indices_of_all_alignments, output, divvy, seen)


'''
This method tells you whether or not a certain divvied result can be found in all other alignments without the need to divvy further. 
If it can be found, you have consensus for that particular divvy, and it can be added to output. 
Otherwise, this divvy does not have full consensus and cannot be written to output. 
'''
def can_this_divvy_be_found_without_subsequent_divvies(toSearch, three_d_array_of_divvies, current_alignment_number):
    for x in range(len(three_d_array_of_divvies)):
        if x == current_alignment_number:
            continue
        divvies_of_current_alignment = three_d_array_of_divvies[x]
        current_alingment_has_match = False
        for divvy in divvies_of_current_alignment:
            match_was_found = True
            for k in range(len(toSearch)):
                current_num = toSearch[k]
                compare_to = divvy[k]
                if current_num == "-":
                    continue
                else:
                    if current_num == compare_to:
                        continue
                    else:
                        match_was_found = False
                        break
            if match_was_found:
                current_alingment_has_match = True
                break
        if current_alingment_has_match == False:
            return False
    return True    
                
    
'''
This function returns a 3d array containing all divvies - the singletons are included by default, but they can be removed easily
'''
def divvies_as_3darray(array_of_dicts):
    return_array = []
    for divvy_dict in array_of_dicts:
        array_per_alignment = []
        for x in divvy_dict:
            #curr_array = divvy_dict[x]
            #number_of_nums = 0
            #for num in curr_array:
            #    if num != "-":
            #        number_of_nums += 1 
            #if number_of_nums >= 2:
            #THE COMMENTED OUT CODE ABOVE IS FOR REMOVING SINGLETONS - HERE WE ARE RETAINING SINGLETONS
            array_per_alignment.append(divvy_dict[x])
        return_array.append(array_per_alignment)
    return return_array


'''
This method returns all divvies as an array of dicts - we will need to further engineer the data returned by this function to make it usable for the purpose of 
writing the correct information to the final output
'''
def get_divvied_results_as_array_of_dicts(all_indices, toSearch):
    someArray = []
    for indices in all_indices:
        myDict = {}
        for i in range(len(indices)):
            if indices[i] != "-":
                myDict[indices[i]] = []
        for i in range(len(indices)):
            if indices[i] == "-":
                for x in myDict:
                    myDict[x].append("-")
            else:
                curr = indices[i]
                myDict[curr].append(toSearch[i])
                for k in myDict:
                    if k != curr:
                        myDict[k].append("-")
        someArray.append(myDict)
    return someArray


'''
In this method, we are passed in a 2D array containing all indices arrays - one indices array per alignment
If none of the indices arrays requires a divvy, return False
If even a single indices array is divvy-able, return True
'''
def divvy_is_required(all_indices):
    for indices_array in all_indices:
        if len(indices_array) == 0:
            continue
        i = 0
        while indices_array[i] == "-":
            i += 1
            continue
        number_to_match_with = indices_array[i]
        for x in range(i, len(indices_array)):
            currentNum = indices_array[x]
            if currentNum == "-":
                continue
            else:
                if currentNum != number_to_match_with:
                     return True
    return False


'''
This method creates the indices array given toSearch, all_alignments, and the current_alignment_number
'''
def create_indices(all_alignments, toSearch, curr_alignment_number):
    indices = []
    current_alignment = all_alignments[curr_alignment_number]
    for x in range(len(current_alignment)):
        row = current_alignment[x]
        num = toSearch[x]
        if num == "-":
            indices.append("-")
        else:
            for a in range(0, len(row)):
                if row[a] == num:
                    indices.append(a)
    return indices


'''
This method creates the first toSearch - this method just considers the first alignment for the purpose of creating a toSearch sequence
'''
def create_toSearch(first_alignment, curr):
    toSearch = []
    for row in first_alignment:
        toSearch.append(row[curr])
    return toSearch


'''
This function stores the transpore of l1 in l2 and returns it
'''
def transpose_output(l1, l2):
    for i in range(len(l1[0])):
        row = []
        for item in l1:
            if i >= len(item):
                break
            row.append(item[i])
        l2.append(row)
    return l2


'''
This function returns the desired output as an arrays of chars - each char is a dash or a letter
The output yielded by this function is processed and written to a text file
'''
def produce_output_in_letters(output, first_alignment):
    return_array = []
    intermediary = []
    for i in range(len(first_alignment[0])):
        string = ""
        for array in first_alignment:
            string += array[i]
        some_array = []
        for char in string:
            if char == "-":
                continue
            some_array.append(char)
        intermediary.append(some_array)
    for nums in output:
        subarray = []
        for i in range(len(output[0])):
            curr_num = nums[i]
            if curr_num == "-":
                subarray.append("-")
            else:
                int_number = int(curr_num) - 1
                subarray.append(intermediary[i][int_number])
        return_array.append(subarray)
    return return_array

'''
This function is passed in a 2d array, and removes each 1d array inside it if it is a singleton 
'''
def remove_singletons(two_d_array):
    return_array = []
    for one_d_array in two_d_array:
        num_letters = 0
        for letter in one_d_array:
            if letter != "-":
                num_letters += 1
        if num_letters > 1:
            return_array.append(one_d_array)
    return return_array


'''
The main function is responsible for data processing (the outputs from the functions provided above need to be engineered a bit) and writing the data to an output file
The output file is formatted as a fasta file
'''
def main():
    [directory_to_use, myDict] = order_the_alignments(alignmentFolder) 
    #directory_to_use contains the correctly ordered fasta files
    #myDict contains the species ids as keys, maps to [species header, species sequence]

    three_d_array_letters = convertToLetters(directory_to_use) #convert fasta files to 3d array of letters

    nums_answer = convert_to_nums(three_d_array_letters) #convert to array of nums and dashes, based on MergeAlign's numbering system
    nums_array_output_initial = create_output(nums_answer)
    nums_array_output = []

    #the for loop below gets rid of duplicates in the number sequence
    for num_array in nums_array_output_initial:
        if num_array not in nums_array_output:
            nums_array_output.append(num_array)
    
    reference_alignment = three_d_array_letters[0] #this is the alignment to which we compare - we recreate the numbering sequence in other alignments based on this array
    letters_array_output_initial = produce_output_in_letters(nums_array_output, transpose_output(reference_alignment, [])) #generates the desired output which needs to be processed
    letters_array_output = letters_array_output_initial
    #letters_array_output = remove_singletons(letters_array_output_initial) #this function removes all singletons before we write our results to output
    letters_array_output = transpose_output(letters_array_output, []) #need to do this for data processing purposes - makes it easier to write to output file later on

    arrayOfSpeciesHeaders = []
    #the for loop below puts all species head in arrayOfSpecies in the correct order
    for key in sorted(myDict.keys()):
        arrayOfSpeciesHeaders.append(myDict[key][0]) #this array contains all species headers in the correct order

    arrayOfSequences = []
    #the for loop below puts all sequences in arrayOfSequences in correct order
    for row in letters_array_output:
        string = ""
        for char in row:
            string += char
        arrayOfSequences.append(string) #this array contains all the sequences in correct order

    with open(last_component+"_result.fasta", "a") as f:
        #this loop writes the output file in fasta format
        for i in range(len(arrayOfSpeciesHeaders)):
            f.writelines(arrayOfSpeciesHeaders[i]) 
            f.writelines(arrayOfSequences[i]+"\n")

main()
