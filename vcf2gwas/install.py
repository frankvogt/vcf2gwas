import subprocess
import os
import shutil

def main():
    #print("Checking environment.. \nUpdating dependencies..")
    #env_name = os.environ["CONDA_DEFAULT_ENV"]
    #subprocess.run(["conda", "env", "update", "-n", env_name, "-f", os.path.join(os.path.dirname(__file__), 'environment.yml')])
    #try:
    #    #subprocess.run(["pip", "install", "umap-learn"])
    #    subprocess.run(["conda", "env", "create", "-f",
    #        os.path.join(os.path.dirname(__file__), 'environment.yml')])
    #except Exception as e:
    #    print(e)
    #print("Environment is now up to date")
    print("Copying example input files to current working directory..")
    try:
        os.mkdir("input")
    except:
        pass
    path = os.path.join(os.path.dirname(__file__), 'input')
    for file in os.listdir(path):
        shutil.copy(os.path.join(path, file), os.path.join('input', file))
    print("Copying README and LICENSE files..")
    shutil.copy(os.path.join(os.path.dirname(__file__), 'README.md'), 'README.md')
    shutil.copy(os.path.join(os.path.dirname(__file__), 'LICENSE'), 'LICENSE')
    #print("Installation complete\n")
