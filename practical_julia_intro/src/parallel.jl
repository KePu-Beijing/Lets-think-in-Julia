using Distributed
using ClusterManagers
using SharedArrays
using ProgressMeter

np = parse(Int, ENV["SLURM_NTASKS_PER_NODE"])
ntask_per_node =  parse(Int, ENV["SLURM_NTASKS_PER_NODE"])
nodes_file = open("nodes_list.log","r")
nodes_list=readlines(nodes_file)
close(nodes_file)

addprocs([("$node-opa",ntask_per_node) for node in nodes_list])

@everywhere using ...

@everywhere begin
  ...
  V_ϕ_tsk =[[i,j] for i in 1:Nϕ for j in i:Nϕ if iseven(i+j)]
  d_V_mtrx=zeros(Nϕ,Nϕ)

  @showprogress pmap(i->V_ϕ(V_ϕ_tsk[i][1],V_ϕ_tsk[i][2], d_V_mtrx), 1:size(V_ϕ_tsk)[1])

  Vmtrx=zeros(Nϕ,Nϕ)
  for i in workers()
      global Vmtrx += remotecall_fetch(()->d_V_mtrx,i)
  end
  ...
end