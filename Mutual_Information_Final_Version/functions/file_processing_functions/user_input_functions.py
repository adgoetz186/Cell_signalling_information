def strip_comments(list_of_lines):
	"""Function that removes comment lines from a list of file lines
			Input:
				list of file lines (list)
			Output:
				list of file lines without comments (list)"""
	list_of_non_comment_lines = []
	for i in list_of_lines:
		if i != "":
			if i[0] != "#":
				list_of_non_comment_lines.append(i)
	return list_of_non_comment_lines

def read_file_to_arg_dict(file_path):
	"""Function that reads a text file specifying user inputs and creates a dictionary of arguments
		Input:
			file path (str)
		Output:
			user specified argument dictionary (dict)"""
	with open(file_path, "r") as user_file:
		user_list = user_file.readlines()
	user_list = strip_comments(user_list)
	user_list = [i.replace(" ", "") for i in user_list]
	return {i.split(":")[0]: eval(i.split(":")[1]) for i in user_list}