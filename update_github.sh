# To generate your own key, see
# https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh
#
# To generate a github key, enter:
# ssh-keygen -t rsa -b 4096 -C "your_github_email@email.com"
# Then you enter the name of the key. Then I add

eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_rsa_github
echo Beginning cloning.
git clone git@github.com:ericgoolsby/GWAS_pipeline.git snuffel_burgers
echo Copying directory.
cp -r snuffel_burgers/* .
echo Deleting directory.
rm -rf snuffel_burgers/
  echo Done.
