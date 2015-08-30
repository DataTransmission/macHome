function tvlct(args)

* SYNTAX  tvlct [ n | 'default' | name ]
*   
*          n - table index
*          'default' - use GrADS rainbow
*          name - table name
*

if ( args = '' )
say 'SYNTAX: tvlct [ n | -default | name | -help ]'
return
endif

if ( args = '-help' )
 say 'SYNTAX: tvlct [ n | -default | name | -help ]'
 say '   -default :  re-load the built-in rainbow sequence'
 say '   n        :  use color table n'
 say '   name     :  use color table by name'
 say '   -help    :  show this message'
 say ' '
 say ' When assigning a color table by name, tvlct will attempt '
 say ' to invoke a procedure by the name  "tvlct_"name.gs'
 say '   i.e. tvlct red-green  will look for tvlct_red-green.gs'
 say ' The pre-defined tables are:'
 say '    1 bigrainbow '
 say '    2 blue_yellow_white'
 say '    3 blue_white_red'
 say '    4 hue_sat_value 1 red-blue-cream'
 say '    5 hue_sat_value 2 periwinkle-green-red'
 say '    6 rainfall colors'
 say ' ' 
endif

if ( args = '-default' | args = 'default' ) 
  'set rbcols'
  return
else; if ( args = '1'  )
    'tvlct_bigrainbow'
    say 'Switching to extended rainbow colors'
    return
else; if ( args = '2'  )
   'tvlct_blue_yellow_red.gs'
   say 'Switching to blue-yellow-red colors'
   return
else; if ( args = '3'  )
   'tvlct_blue_white_red.gs'
   say 'Switching to blue-white-red colors'
   return
else; if ( args = '4'  )
   'tvlct_hsv1.gs'
   say 'Switching to hue_sat_value1 colors'
   return
else; if ( args = '5'  )
   'tvlct_hsv2.gs'
   say 'Switching to hue_sat_value2 colors'
   return
else; if ( args = '6'  )
   'tvlct_rainfall.gs'
   say 'Switching to rainfall colors'
   return
else
    'tvlct_'args'.gs'
    return
endif
endif
endif
endif
endif
endif
endif

cmd='tvlct_'args



  
