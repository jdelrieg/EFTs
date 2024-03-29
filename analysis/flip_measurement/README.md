# flip_measurement 

This directory contains several scripts that are relevant for the charge flip measurement.

### Processors 

* `flip_mr_processor.py`: This processor is responsible for performing the measurement of the flip probabilities. It is designed to process MC, and based on the MC truth info, counts the number of electrons whose charge has been mismeasured. The information is saved into a 2d `coffea` histogram (with axes for `pt` and `abs(eta)`). 

* `flip_ar_processor.py`: This processor can be used for validations purposes. It expects to process data. It applies the measured flip probabilities to OS events in order to estimate the number of SS events. The region it uses for this application test is the flip CR that we have added to `topeft`. 

### Plotters

* `flip_mr_plotter.py`: This script is designed to run on the output histograms of `flip_mr_processor.py`. It plots the 2d pt-abs(eta) histograms, and also saves these histograms into pkl files. The pkl files can then be copied into `topcoffea/data/fliprates` to be used by `corrections.py`.

* `flip_ar_plotter.py`: This script is designed to plot the output histograms of `flip_ar_processor.py`. It is useful as a simple comparison between the SS data and the prediction (though using `topeft` and looking at the flip CR would be a much more thorough comparison since that would incorporate other contributions into the prediction, e.g. `fakes`). 

### Run scripts

* `run_flip.py`: This run script can run either the `flip_mr_processor.py` or `flip_mr_processor.py` processor (specify the name of the processor with the `--processor_name` argument). It can also run with the `futures` or `work_queue` executor (specify this via the `--executor` argument).
