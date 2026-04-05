from TB2J_OpenMX.ffiparser import OpenmxWrapper
from TB2J.exchange import ExchangeNCL, ExchangeCL
import os


def gen_exchange(path='./', prefix='openmx', **kwargs):
    tbmodel = OpenmxWrapper(path, prefix)

    # efermi is read from OpenMX data, not user input
    kwargs['efermi'] = tbmodel.efermi

    # Auto-generate description and prepend to any user-provided one
    auto_desc = f"""Using OpenMX data: 
path: {os.path.abspath(path)}
prefix: {prefix}
"""
    user_desc = kwargs.get('description') or ''
    kwargs['description'] = auto_desc + user_desc

    ExchangeClass = ExchangeNCL if tbmodel.non_collinear else ExchangeCL

    print("Starting to calculate exchange.")
    exchange = ExchangeClass(
        tbmodels=tbmodel,
        atoms=tbmodel.atoms,
        basis=tbmodel.basis,
        **kwargs,
    )
    output_path = kwargs.get('output_path', 'TB2J_results')
    exchange.run(path=output_path)
    print(f"\nAll calculation finished. The results are in {output_path} directory.")

if __name__=='__main__':
    gen_exchange(
        path='/home/hexu/projects/TB2J_example/OPENMX/SrMnO3_FM_SOC/', magnetic_elements=['Mn'], nz=50, Rcut=8)
