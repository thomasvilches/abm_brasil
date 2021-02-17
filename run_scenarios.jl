 # states = ["saopaulo";"riograndedosul";"bahia";"para";"goias"]
states = ["saopaulo"]

for city in states
   #=  run_param_fixed(city,[10;20;30],1.0,1.0)

    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,0.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,0.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,0.0)

    ###supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,1.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,1.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,1.0) 
 

    #### All the same scenarios, supposing that who got vaccinated does not isolate themselve


    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,1.0,0.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,1.0,0.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,1.0,0.0)

    ###supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,1.0,1.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,1.0,1.0)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,1.0,1.0)

    ##DOUBLEVAC
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,0.0,600)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,0.0,600)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,0.0,600)

    ###supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,1.0,600)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,1.0,600)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,1.0,600)  

 =#
  
    ####OXFORD
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,0.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,0.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,0.0,300,1,0.641,84,[21;14])

    ###supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,1.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,1.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,1.0,300,1,0.641,84,[21;14]) 
 

    ###OXFORD
    #### All the same scenarios, supposing that who got vaccinated does not isolate themselve


    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,1.0,0.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,1.0,0.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,1.0,0.0,300,1,0.641,84,[21;14])

    ###supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,1.0,1.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,1.0,1.0,300,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,1.0,1.0,300,1,0.641,84,[21;14])
 
   ###DOUBLEVAC
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,0.0,600,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,0.0,600,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,0.0,600,1,0.641,84,[21;14])

    ###supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,1.0,600,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,1.0,600,1,0.641,84,[21;14])
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,1.0,600,1,0.641,84,[21;14])

#= 
    #### scenario vac_ef = 50/70%
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,1.0,0.0,0.0,300,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.5,0.0,0.0,300,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.0,0.0,0.0,300,1,0.35,28)

    ###supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,1.0,0.0,1.0,300,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.5,0.0,1.0,300,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.0,0.0,1.0,300,1,0.35,28) 
 

    #### All the same scenarios, supposing that who got vaccinated does not isolate themselve


    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,1.0,1.0,0.0,300,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.5,1.0,0.0,300,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.0,1.0,0.0,300,1,0.35,28)

    ###supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,1.0,1.0,1.0,300,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.5,1.0,1.0,300,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.0,1.0,1.0,300,1,0.35,28)

    ##DOUBLEVAC
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,1.0,0.0,0.0,600,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.5,0.0,0.0,600,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.0,0.0,0.0,600,1,0.35,28)

    ###supposing that vaccine protects 100% agains severe disease
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,1.0,0.0,1.0,600,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.5,0.0,1.0,600,1,0.35,28)
    run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.0,0.0,1.0,600,1,0.35,28)  =#
end 


 

a = 5:10:45
b = [i for i in a] 
#b = [45]
for city in states
    for i in b
         run_param_fixed(city,[10;20;30],1.0,1.0,false,0.0,0.0,0.0,0.0,300,i)

        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,0.0,300,i)
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,0.0,300,i)
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,0.0,300,i)

        ###supposing that vaccine protects 100% agains severe disease
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,1.0,300,i)
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,1.0,300,i)
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,1.0,300,i) 
    

        ####OXFORD
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,0.0,300,i,0.641,84,[21;14])
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,0.0,300,i,0.641,84,[21;14])
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,0.0,300,i,0.641,84,[21;14])

        ###supposing that vaccine protects 100% agains severe disease
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,1.0,300,i,0.641,84,[21;14])
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,1.0,300,i,0.641,84,[21;14])
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,1.0,300,i,0.641,84,[21;14])   

        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,1.0,0.0,0.0,300,i,0.35,28)
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.5,0.0,0.0,300,i,0.35,28)
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.0,0.0,0.0,300,i,0.35,28)
    
        ###supposing that vaccine protects 100% agains severe disease
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,1.0,0.0,1.0,300,i,0.35,28)
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.5,0.0,1.0,300,i,0.35,28)
        run_param_fixed(city,[10;20;30],1.0,1.0,true,0.7,0.0,0.0,1.0,300,i,0.35,28) 
     
    
    
    end
end