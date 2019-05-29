import subprocess


# def alreadyExistFolderName(folderNameToSave, path):
#     # filenameToSave become unique
#     for root, dirs, files in os.walk(path):
#         for folder in dirs:
#             folder = folder + "/"
#             if(folder == folderNameToSave):
#                 return True
#     return False


# # find path and complete name of the folder
# def findOutputPathFolder(path):
#     currentAttempt = 1
#     titleAttempt = '{0:03}'.format(currentAttempt)
#     folderNameToSave = str(currentDate) + "_" + str(titleAttempt) + "/"

#     # if name exist increment titleAttempt
#     while(alreadyExistNameFile(folderNameToSave, path)):
#         currentAttempt += 1
#         currentTitleAttempt = '{0:03}'.format(currentAttempt)
#         folderNameToSave = str(currentDate) + "_" + \
#             str(currentTitleAttempt) + "/"

#     outputFolderPathName = path + folderNameToSave

#     return outputFolderPathName


if __name__ == "__main__":
    # path = "../result/R05_hmmerTest"
    # outputPathFolderName = findOutputPathFolder(path)
    # if not os.path.exists(outputPathFolder):
    #     os.makedirs(outputPathFolder)

    # PER FUNZIONARE BISOGNA AVERE LA CARTELLA BiologicalDataHmmer CON DENTRO I FILES DI HMMER

    # subprocess.call(['cd..'], shell=True)
    # subprocess.call(['cd ./BiologicalDataHmmer/hmmer-3.2.1'], shell=True)
    # subprocess.call(['bash'], shell=True)
    # subprocess.call(['dir'], shell=True)

    # subprocess.call(
    #     ['hmmbuild {} {}'.format(actualOutputPathName,actualInputMsaCleaned)], shell=True)

    def bash_command(cmd):
        subprocess.Popen(cmd, shell=True, executable='bash')

    bash_command('ls')
