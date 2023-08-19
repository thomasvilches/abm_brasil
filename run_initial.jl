states = ["saopaulo"]
a = 5:10:45
b = [i for i in a] 
#b = [45]
for city in states
    for i in b
        run_param_initial(city,[10;20;30],1.0,1.0,false,0.0,0.0,0.0,0.0,300,i)

        run_param_initial(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,0.0,300,i)
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,0.0,300,i)
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,0.0,300,i)

        ###supposing that vaccine protects 100% agains severe disease
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.5038,1.0,0.0,1.0,300,i)
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.5038,0.5,0.0,1.0,300,i)
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.5038,0.0,0.0,1.0,300,i) 
    

        ####OXFORD
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,0.0,300,i,0.641,84,[21;14])
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,0.0,300,i,0.641,84,[21;14])
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,0.0,300,i,0.641,84,[21;14])

        ###supposing that vaccine protects 100% agains severe disease
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.7042,1.0,0.0,1.0,300,i,0.641,84,[21;14])
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.7042,0.5,0.0,1.0,300,i,0.641,84,[21;14])
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.7042,0.0,0.0,1.0,300,i,0.641,84,[21;14])  
        
            
        #*###PFIZER
        #! No protection severe
        #https://docs.google.com/document/d/1sIObk1Uaq0_ja2spTkeus1hO-HNjhiA63d4owty-eQs/edit
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.94,1.0,0.0,0.0,300,i,0.57,60,[14;7])
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.94,0.5,0.0,0.0,300,i,0.57,60,[14;7])
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.94,0.0,0.0,0.0,300,i,0.57,60,[14;7])

        #!##supposing that vaccine protects 100% agains severe disease
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.94,1.0,0.0,1.0,300,i,0.57,60,[14;7])
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.94,0.5,0.0,1.0,300,i,0.57,60,[14;7])
        run_param_initial(city,[10;20;30],1.0,1.0,true,0.94,0.0,0.0,1.0,300,i,0.57,60,[14;7]) 

        # run_param_initial(city,[10;20;30],1.0,1.0,true,0.7,1.0,0.0,0.0,300,i,0.35,28)
        # run_param_initial(city,[10;20;30],1.0,1.0,true,0.7,0.5,0.0,0.0,300,i,0.35,28)
        # run_param_initial(city,[10;20;30],1.0,1.0,true,0.7,0.0,0.0,0.0,300,i,0.35,28)
    
        # ###supposing that vaccine protects 100% agains severe disease
        # run_param_initial(city,[10;20;30],1.0,1.0,true,0.7,1.0,0.0,1.0,300,i,0.35,28)
        # run_param_initial(city,[10;20;30],1.0,1.0,true,0.7,0.5,0.0,1.0,300,i,0.35,28)
        # run_param_initial(city,[10;20;30],1.0,1.0,true,0.7,0.0,0.0,1.0,300,i,0.35,28) 
     
    
    
    end
end