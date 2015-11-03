#export PATH=/usr/local/cuda/bin:$PATH
export PATH=/usr/local/bin:$PATH
export PATH=/Users/ginochen/scripts/shscript:$PATH
export PATH=~/Documents/software/netcdf-4.1.3/ncdump/:$PATH
export PATH=/usr/texbin/:$PATH
export PATH=/usr/bin/:$PATH
#export PATH=/Users/ginochen/Downloads/pdfsizeopt/:$PATH
export DYLD_LIBRARY_PATH=/usr/local/cuda/lib:$DYLD_LIBRARY_PATH
export GADDIR=~/software/grads-1.9b4/data/
export EDITOR=vim
# Make ls colorful in the following two lines
export CLICOLOR=1
export LSCOLORS=ExFxCxDxBxegedabagacad

# For Go using Google Drive
export GOROOT=/usr/local/go
export GOPATH=$HOME/gopath
export PATH=$GOPATH:$GOPATH/bin:$GOROOT:$GOROOT/bin:$PATH


alias lsf='ls -F'
alias safari='open /Applications/Safari.app/'
alias iphone='ssh -XY root@192.168.1.69'
alias mat='/Applications/MATLAB_R2015a.app/bin/matlab -nodesktop -nosplash'
alias mkdirp='mkdir -p' # make a full tree path ($ mkdirp /this/is/a/tree/path)
alias tarf='tar -cvf' # tar -cvf filename.tar files*
alias grep='grep --color=always -n -r'
alias preview='open -a Preview'
alias webmount='mount -t smbfs //gchen@krystal.rsmas.miami.edu/gchen ~/WebShare'
alias webumount='umount ~/WebShare'
alias tree='tree -C'

homedir='/Users/ginochen/'

#export CDPATH=$homedir/notes # use cautiously since this will have conflict when path name not under CDPATH has the same name

alias cdqj='cd $homedir/Gino/Paper/Paper_QJRMS/paper_revision2/submitted/Supplementary_Material_not_for_Review/Tex_documents/'
alias cdnote='cd $homedir/notes'
alias cdconvex='cd $homedir/notes/convex'
alias cdlinearalg='cd $homedir/notes/linearAlg'
alias cdbook='cd $homedir/books'
alias cdpaper='cd $homedir/paper'
alias cdmvim='cd $homedir/vimfiles/MacVim-snapshot-74/./mvim'
alias cdreciept='cd $homedir/Documents/Reciept'

#alias ..="cd .."
#alias ..2="cd ../.."
#alias ..3="cd ../../.."
#alias ..4="cd ../../../.."
#alias ..5="cd ../../../../.."
alias cd..="cd .." # navigate up 1
alias cd...="cd ../.." #navigate up 2, etc
alias cd....="cd ../../.."
alias cd.....="cd ../../../.."
alias cd......="cd ../../../../.."

# ssh shortcuts
# -C  Requests compression of all data (including stdin, stdout, stderr, and data for forwarded X11 and TCP/IP connections).
# -Y  Enables trusted X11 forwarding.
alias sshenso='ssh -Y -C gchen@enso.rsmas.miami.edu' 
alias sshsnowy='ssh -Y -C gchen@snowy.rsmas.miami.edu'
alias sshares='ssh -Y -C gchen@ares.ccs.miami.edu'
alias sshkaos='ssh -Y -C gchen@kaos.rsmas.miami.edu'
alias sshpeg2='ssh -Y -C gchen@pegasus2.ccs.miami.edu'
alias sshvisx='ssh -Y -C gchen@visx.ccs.miami.edu' #10.141.249.4
alias sshys='ssh -YC -l ginochen yellowstone.ucar.edu'
alias sshmko='ssh -Y -C mko@129.171.98.195'


#------------------------
# Github section
#------------------------
# CREATE LOCAL REPO: a new local repo named "gino" and sync to the remote https://github.com/ginochen/gino.git with the short name "origin"
# git remote add origin https://github.com/ginochen/gino.git
# PUSH TO REMOTE REPO: push the local branch "master" to the shortname "origin" on the remote repo
# git push origin master
# BACK TO OLD VERSION:  show old versions and checkout to old versions
# git log; git checkout versionnumber
alias gitcommit='git commit -a' # automatically stage the local files, don't need to go to editor 
# 1. Goto local branch ~/Github/gino
# 2. Commit every changes without adding new files
# 3. Push everything into the remote repo https://github.com/ginochen/gino.git under the master folder ~/Github/gino
#    If the master folder is ~/Github/Hsiaochu, then the corresponding remote repo becomes https://github.com/ginochen/Hsiaochu.git
alias gitpush='git push -u origin master'
alias gitpush2pegMaster='git push ssh://gchen@pegasus2.ccs.miami.edu:/nethome/gchen/macHome.git master'
alias gitpush2pegBranch='git push ssh://gchen@pegasus2.ccs.miami.edu:/nethome/gchen/macHome.git macHome2peg' # this branch tracks large data (e.g., books)
alias gitrm='git rm -cached'

alias untar='tar xvfz'

# there is no more ~/Github folder in macHome dir so change this line
alias rsyncTouro='rsync -av ~/Desktop/ /Volumes/Hsiaochu_Gino/Desktop; rsync -av ~/Github/ /Volumes/Hsiaochu_Gino/Github; rsync -av ~/Documents/ /Volumes/Hsiaochu_Gino/Documents'



# SRINK pdf size with jbig2 as a pointer to the same words
alias pdfsizeopt='~/Downloads/pdfsizeopt/./pdfsizeopt --use-pngout=false --use-jbig2=true --use-multivalent=false' 
# or transfer .djvu file into .tif with '$ddjvu -format=tif -quality=75 file.djvu filename.tif' and use '$ jbig2 -b filename -p -s filename.tif', 
# lastly use '$ pdf.py filename > filename.pdf' to turn into pdf file

# shortcut for directories
function tabname {
  printf "\e]1;$1\a"
}
 
function winname {
  printf "\e]2;$1\a"
}
##
# Your previous /Users/ginochen/.profile file was backed up as /Users/ginochen/.profile.macports-saved_2012-06-05_at_16:54:26
##


# Useful functions
# make a long path and cd to it at the same time
function mkdircd () { mkdir -p "$@" && eval cd "\"\$$#\""; } # mkdir -p /this/is/a/long/path; cd /this/is/a/long/path




# MacPorts Installer addition on 2012-06-05_at_16:54:26: adding an appropriate PATH variable for use with MacPorts.
#export PATH=/opt/local/bin:/opt/local/sbin:$PATH
# Finished adapting your PATH environment variable for use with MacPorts.


##
# Your previous /Users/ginochen/.profile file was backed up as /Users/ginochen/.profile.macports-saved_2014-11-15_at_18:20:29
##

# MacPorts Installer addition on 2014-11-15_at_18:20:29: adding an appropriate PATH variable for use with MacPorts.
export PATH="/opt/local/bin:/opt/local/sbin:$PATH"
# Finished adapting your PATH environment variable for use with MacPorts.


##
# Your previous /Users/ginochen/.profile file was backed up as /Users/ginochen/.profile.macports-saved_2015-06-06_at_09:28:41
##

# MacPorts Installer addition on 2015-06-06_at_09:28:41: adding an appropriate PATH variable for use with MacPorts.
export PATH="/opt/local/bin:/opt/local/sbin:$PATH"
# Finished adapting your PATH environment variable for use with MacPorts.

