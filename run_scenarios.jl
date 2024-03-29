 # states = ["saopaulo";"riograndedosul";"bahia";"para";"goias"]
states = ["saopaulo"]

for city in states
    #* No vac escenario
    run_param_fixed(city,[10;20;30],1.0,1.0)

    #* CORONAVAC
    #! No protection severe
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,0.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,0.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,0.0)

    #! ##supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,1.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,1.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,1.0) 
 

    #?### All the same scenarios, supposing that who got vaccinated does not isolate themselve

    #! No protection severe
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,1.0,0.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,1.0,0.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,1.0,0.0)

    #!##supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,1.0,1.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,1.0,1.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,1.0,1.0)

    #?#DOUBLEVAC
    #! No protection severe
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,0.0,600)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,0.0,600)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,0.0,600)

    #!##supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,1.0,600)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,1.0,600)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,1.0,600)  

 
  
    #*###OXFORD
    #! No protection severe
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,0.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,0.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,0.0,300,1,0.641,84,[21;14])

    #!##supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,1.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,1.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,1.0,300,1,0.641,84,[21;14]) 
 

    #?### All the same scenarios, supposing that who got vaccinated does not isolate themselve
    #! No protection severe
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,1.0,0.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,1.0,0.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,1.0,0.0,300,1,0.641,84,[21;14])

    #!##supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,1.0,1.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,1.0,1.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,1.0,1.0,300,1,0.641,84,[21;14])
 
    #?##DOUBLEVAC
    #! No protection severe
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,0.0,600,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,0.0,600,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,0.0,600,1,0.641,84,[21;14])

    #!##supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,1.0,600,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,1.0,600,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,1.0,600,1,0.641,84,[21;14])

    
    #*###PFIZER
    #! No protection severe
    #https://docs.google.com/document/d/1sIObk1Uaq0_ja2spTkeus1hO-HNjhiA63d4owty-eQs/edit
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,1.0,0.0,0.0,300,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.5,0.0,0.0,300,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.0,0.0,0.0,300,1,0.57,60,[14;7])

    #!##supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,1.0,0.0,1.0,300,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.5,0.0,1.0,300,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.0,0.0,1.0,300,1,0.57,60,[14;7]) 
 

    #?### All the same scenarios, supposing that who got vaccinated does not isolate themselve
    #! No protection severe
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,1.0,1.0,0.0,300,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.5,1.0,0.0,300,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.0,1.0,0.0,300,1,0.57,60,[14;7])

    #!##supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,1.0,1.0,1.0,300,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.5,1.0,1.0,300,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.0,1.0,1.0,300,1,0.57,60,[14;7])
 
    #?##DOUBLEVAC
    #! No protection severe
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,1.0,0.0,0.0,600,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.5,0.0,0.0,600,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.0,0.0,0.0,600,1,0.57,60,[14;7])

    #!##supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,1.0,0.0,1.0,600,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.5,0.0,1.0,600,1,0.57,60,[14;7])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.94,0.0,0.0,1.0,600,1,0.57,60,[14;7])

end 