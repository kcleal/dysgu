from subprocess import run
import os
import click


@click.command()
@click.option('-j', default=1, help='Number of threads')
def build(j):

    print("Building htslib")
    sub = "./dysgu/htslib"
    if not os.path.exists(sub):
        print(f"htslib is missing from {sub}, try downloading with git clone --recursive https://github.com/kcleal/dysgu.git")
        quit()
    run(f"cd ./dysgu/htslib; autoheader; autoconf; ./configure; make -j {j};", shell=True)

    print("Installing dependencies")
    run(["pip install -r requirements.txt"], shell=True)

    print("Installing dysgu")
    run(["python setup.py install"], shell=True)

    run(["dysgu --version"], shell=True)
    run(["dysgu test"], shell=True)
    print("Done")


if __name__ == '__main__':
    build()
