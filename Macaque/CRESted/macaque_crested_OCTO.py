import enhancer_suite_crested as esc

# def main():

# load example data configuration
crested_data_config = esc.CREstedDataConfig(
    data_root = "/home/nelson/bg_macaque_model_data",
    dataset = "basalganglia",
    project_name = "macaque_basalganglia_crested_filtered",
    run_name = "crested_bg_filtered",
    species="macaque",
    sequence_length=2114,
    training_epochs=100,
    batch_size=64,
)
print(f"CREsted Data Configuration: {crested_data_config}")

# load example model configuration
crested_model_config = esc.CREstedModelConfig(
    model_name="ChromBPNet",
    sequence_length=2114,
    num_classes=48,
)
print(f"CREsted Model Configuration: {crested_model_config}")

# instantiate example model
crested_model = esc.CREstedModel(config=crested_model_config)

# execute model training
crested_model.train(crested_data_config)


if __name__ == "__main__":
    main()