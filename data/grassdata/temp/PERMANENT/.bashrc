test -r ~/.alias && . ~/.alias
PS1='GRASS 6.4.2 (temp):\w > '
PROMPT_COMMAND="'/usr/lib64/grass-6.4.2/etc/prompt.sh'"
export PATH="/usr/lib64/grass-6.4.2/bin:/usr/lib64/grass-6.4.2/scripts:/home/erodrigu/.grass6/addons:/home/erodrigu/anaconda/envs/SWOTRiver/bin:/home/erodrigu/sw/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/home/erodrigu/bin"
export HOME="/home/erodrigu"
export GRASS_SHELL_PID=$$
trap "echo \"GUI issued an exit\"; exit" SIGQUIT
