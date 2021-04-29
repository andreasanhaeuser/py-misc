# py-misc
Various python utils


How to install
==============
1) Go to a folder where you want to repo to be installed. On Andi's computer,
this is ~/shared/lib but this can be any folder on your system.
1) Make sure this folder is in your $PYTHONPATH
1) Make sure you have access to the repo (it's private, ask Andi or Aidan if
you haven't)
1) Create an ssh keypair
1) Upload your public key to your github profile:
    * On the github web interace click on your avatar (top right)
    * Go to "Settings"
    * Go to SSH and GPG keys
    * Click "New SSH key" and upload your PUBLIC (!) key
1) It may be necessary to specify your credentials in your ~/.ssh/config
   (change the `IdentityFile` entry to match the path to you ssh key):
    Host github.com
      HostName github.com
      User git
      IdentityFile /home/andreas/.ssh/github
      IdentitiesOnly yes
      PreferredAuthentications publickey
1) Clone the repository to this folder under the name `misc`
    (note the missing `py-`):
   git clone git@github.com:/andreasanhaeuser/py-misc misc
