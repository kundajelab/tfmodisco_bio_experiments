hyperas_opt_gpu1.py   100 trials, distributed using mongodb localhost:12345, running in ~/mongodb, ipynb screen
                      mongod --dbpath . --port 12345 &   

master: surya:gpu1,   
                      $TFNET_ROOT/scripts/run_pipeline.py --tfs SPI1 --end-task 4 --start 40 --end 50
                      hyperas_prepare.py
                      python ../hyperas_opt_gpu1.py &> opt1.txt"

                      mongoexport --host surya --port 12345 --db foo_db -c jobs --csv -f result.model_param  -o db.csv

worker:               copy over temp_model.py, data_big.npz
                      hyperopt-mongo-worker --mongo=localhost:12345/foo_db --exp-key='exp1' &> worker1.txt

workers:              surya:gpu1, gpu2, 
                      nandi:gpu1(CUDA3), tfnet(CUDA2), gpu5(night), gpu7(night), 
                      kat-pod1a:gpu0, gpu1, 
                      kat-pod1b:gpu0,gpu1


