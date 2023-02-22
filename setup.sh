pathJulia=""
# needs to be adjusted
customPathJulia=/home/stephan/Code/julia/julia
versionJulia=""
requiredJuliaVersion="1.8.5"

# first try system julia
if [ ! -z "$(whereis julia | cut -d ':' -f 2)" ]; then
   echo "Using system's julia if available"
   pathJulia=$(whereis julia | cut -d ':' -f 1)
   versionJulia=$(whereis julia | cut -d ':' -f 2)
fi

# if version not the required, use custom julia installation
if [ "$versionJulia" != "$requiredJuliaVersion" ]; then
   echo -n "System's Julia outdated, trying custom julia installation..."
   if [ ! -z "$customPathJulia" ]; then
      pathJulia=$customPathJulia
      versionJulia=$($pathJulia --version | awk -F ' ' '{print $3}')
   else
      echo "No custom Julia provided. Exiting."
   fi
else
   echo "System's Julia... "
fi

if [ "$versionJulia" != "$requiredJuliaVersion" ]; then
   echo "Neither system nor custom julia installation version satisfied. Exiting."
   exit 1
else
   echo " found!"
fi

# Set up Julia PEtab package
git clone https://github.com/CleonII/Master-Thesis.git
cd Master-Thesis
echo "Installing needed Python packages"
echo "Installing Julia packages (might take some time https://xkcd.com/303/)"
${pathJulia} --project=. -e"using Pkg; Pkg.instantiate()"
echo "Done"

 Set up PyPesto and AMICI
cd ../pypesto_benchmark
bash setup.sh

cd ..
echo "Done setting up"
