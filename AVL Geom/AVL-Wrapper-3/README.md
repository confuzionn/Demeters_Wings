# AVL-Wrapper
### AVL wrapper
ADaX stands for AVL Designer and Executor. It is a Python library and mini-interface for creating ".avl" and ".run" files and automatically running them through AVL.

### Structure
There are a few files of interest, all located in the top-level folder. In descending order, they are:

- "ADaX.py" is the main script to run when you want to run AVL.

- "Aircraft_Geometry.py" contains information related to creating/updating .avl files. It may be customized to create your plane geometry by
  defining everything inside the Update_Geometry() function or by accepting keyword arguments (kwargs) to loop through geometry configurations.

- "AVL_Presets.py" is a file containing some commonly used run cases. Some presets must be modified to fit your specific .avl file, 

- "avl.exe" is the actual AVL application, written in Fortran and compiled into a .exe. You may open the executable and play with AVL in the 
  command line.

- "CaseMaker.py" may be used to create custom run cases for easily storing and running a run case. This is not necessary to set up until later
  in the design process; it may be desired for testing many configurations for many flight phases.

- "Log.txt" is AVL's command line output written into a .txt file, mainly for debugging, but outputs may be read from Log.txt.

- "Raw Run Commands" is a extension-less file where the raw AVL commands may be put into a text file. No fancy commands from the ADaX library,
  just a simple AVL wrapper.

### Documentation
The wrapper code is stored in "Code/Codebase.py", and most functions have a short description. Documentation is in the Documentation folder, including an AVL introduction, examples of real ADaX use, cheatsheets, and guides to the files listed above.

### Getting Started
Open AVL-Wrapper-main in Visual Studio, as all the code's directories are based on this folder being the working directory. Maybe someone in the future could do some directory management. 

### Future Work and Developers
There is room for improvement in many areas of this wrapper. Some improvements are listed below:
 - A GUI interface through the tkinter library (very big undertaking)
 - Directory management
 - Memory management (not sure if it'll ever matter)
 - More descriptive variable names
 - Run time improvements to sweep()
 - Removing the developers' unhinged 3 a.m. comments

If you are interested in making permanent changes to the wrapper, make sure to update documentation and cheatsheets.
Feel free to contact the developers or testers for questions on the wrapper. 

# Current Developers/Testers:
| **Name** | **Discord (AMAT Nickname)** | 
|-|-|
|Cameron M| stowon (Cameron M)|
|Esther K| estherberry (Esther)|



