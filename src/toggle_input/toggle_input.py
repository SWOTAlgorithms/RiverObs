
"""
Toggles the input cells of a python notebook.

To use, from the notebook, execute this cell:

from IPython.display import HTML
from toggle_input import toggle_input
HTML(toggle_input)
"""
##from IPython.display import HTML

##HTML('''<script>
##code_show=true; 
##function code_toggle() {
## if (code_show){
## $('div.input').hide();
## } else {
## $('div.input').show();
## }
## code_show = !code_show
##} 
##$( document ).ready(code_toggle);
##</script>
##The raw code for this IPython notebook is by default hidden for easier reading.
##To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>.''')
##

toggle_input = '''<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
The raw code for this IPython notebook is by default hidden for easier reading.
To toggle on/off the raw code, click <a href="javascript:code_toggle()">here</a>.'''

