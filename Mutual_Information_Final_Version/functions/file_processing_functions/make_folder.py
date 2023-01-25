import os
import shutil


def make_folder(folder_dir,ask_user = True):
	"""Creates folder at the specified location, if folder of same name already exists, gives option for deletion
	returns an error if folder contains sub-folders

		Input:

			-folder_dir (str): location and name of folder to create

			-ask_user (boolean): if True takes effect if there is a folder which already exists at the specified
			location. Prompts the user to specify if they want the pre existing file to be deleted.

		Output:

			-returns none but generates a folder at the specified location"""
	try:
		# tries to make folder
		os.mkdir(folder_dir)
	except OSError as error:
		# excepts OS Error
		# Checks to see if folder contains sub folders
		cont_dir = False
		for i in os.listdir(folder_dir):
			if os.path.isdir(folder_dir+"/"+i):
				cont_dir = True
		# Triggers if folder contains sub folders. This func. should not be used on directories with sub-directories.
		# This is a precaution to help make sure passing the wrong string doesn't delete important directories
		if cont_dir:
			raise OSError
		# Prompts user for input
		if ask_user:
			delete_file = input(
				f"{folder_dir} already exists, remove? (y/n)")
		else:
			delete_file = 'y'
		# if user gives permission, fully deletes folder and creates a new one
		if delete_file.lower() == "y":
			shutil.rmtree(f"{folder_dir}")
			os.mkdir(f"{folder_dir}")
		else:
			raise OSError