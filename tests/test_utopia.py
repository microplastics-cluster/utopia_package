from utopia.utopia import utopiaModel
from utopia.results_processing.mass_balance_check import *
from utopia.results_processing.process_results import *

# Create the model, pass it example config and data, then run it
config_data = utopiaModel.load_json_file("data/default_config.json")
data_data = utopiaModel.load_json_file("data/default_data.json")
model = utopiaModel(config=config_data, data=data_data)
model.summarize()
model.run()

# Process results and print and plot them
massBalance(model)
processor = ResultsProcessor(model)
processor.estimate_flows()
processor.generate_flows_dict()
processor.process_results()
for fraction in ["mass_fraction", "number_fraction"]:
    processor.plot_fractionDistribution_heatmaps(fraction)

processor.extract_results_by_compartment()
for fraction in ["%_mass", "%_number"]:
    processor.plot_compartment_distribution(fraction)
processor.results_by_comp

for i in range(len(processor.results_by_comp)):
    emissions = sum(
        processor.model.emiss_dict_g_s[
            processor.results_by_comp["Compartments"].iloc[i]
        ].values()
    )
    print(
        f"Mass balance for {processor.results_by_comp['Compartments'].iloc[i]}: {processor.results_by_comp['Total_inflows_g_s'].iloc[i]+emissions-processor.results_by_comp['Total_outflows_g_s'].iloc[i]}"
    )

# Calculate exposure indicators
processor.estimate_exposure_indicators()
processor.overall_exposure_indicators
processor.size_fraction_indicators
processor.estimate_emission_fractions()
