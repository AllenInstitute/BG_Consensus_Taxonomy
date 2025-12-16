import enhancer_suite_crested as esc

# def main():

# load example data configuration
crested_data_config = esc.CREstedDataConfig(
    data_root = "/home/nelson/bg_human_model_data",
    dataset = "basalganglia",
    project_name = "human_basalganglia_crested_filtered",
    run_name = "crested_bg_filtered",
    species="human",
    sequence_length=2114,
    training_epochs=100,
    batch_size=64,
)
print(f"CREsted Data Configuration: {crested_data_config}")

# load example model configuration
crested_model_config = esc.CREstedModelConfig(
    model_name="ChromBPNet",
    sequence_length=2114,
    num_classes=52,
)
print(f"CREsted Model Configuration: {crested_model_config}")

# instantiate example model
crested_model = esc.CREstedModel(config=crested_model_config)

# execute model training
crested_model.train(crested_data_config)


if __name__ == "__main__":
    main()



# beaker session create \
#   --budget ai2/ai1 \
#   --gpus 4 \
#   --workspace ai2/enhancer_suite \
#   --image beaker://kasiak/enhancer_suite \
#   --mount src=hostpath,ref=/home/nelson/github/enhancer_suite,dst=/code/enhancer_suite