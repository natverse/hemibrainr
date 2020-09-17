# read hemibrain meshes

hemibrain_read_neuronmesh_single <-function(x){

  neuron_mesh = read_cloudvolume_meshes(segments = x,
                                        cloudvolume.url='precomputed://gs://neuroglancer-janelia-flyem-hemibrain/segmentation_52a13')

}
