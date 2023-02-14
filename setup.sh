
# To stay consistent let us use Julia 1.8.5
pathJulia="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia"
if [ ! "$(eval "${pathJulia} --version")"  == "julia version 1.8.5" ];then
    >&2 echo "Error : Julia path most point to Julia 1.8.5"
    exit 1
fi

# To use Conda (Fides in Julia)
eval "$(conda shell.bash hook)"

# Set up Julia PEtab package
git clone https://github.com/CleonII/Master-Thesis.git
cd Master-Thesis
echo "Installing needed Python packages"
conda env create -f PeTab.yml
echo "Done"
conda activate PeTab
echo "Installing Julia packages (might take some time https://xkcd.com/303/)"
${pathJulia} --project=. -e"using Pkg; Pkg.instantiate()"
echo "Done"

# Set up PyPesto and AMICI
cd ../pypesto_benchmark
bash setup.sh

cd ..
echo "Done setting up"
