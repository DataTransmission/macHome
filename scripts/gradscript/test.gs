'open obsst_201008_plev.ctl'
'open sys3_201008_plev.ctl'
'set gxout shaded'
count = 20
tn = 61
while (count < tn)
   'set t 'count
   prompt 'press enter to plot day ' math_nint(count/2)
*round up tp closest integer
   pull s
   'subplot 2 1 1'
   'd q.1'
   'subplot 2 2 1'
   'd q.2'
   'printim out' count 'jpg' 
   count = count + 1
endwhile 

