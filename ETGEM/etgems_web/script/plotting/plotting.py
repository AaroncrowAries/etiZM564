import matplotlib.pyplot as plt
import numpy as np

#Visualization of calculation results
def Draw_Biomass_By_Glucose_rate(Biomass_list,xmin,xmax,step,ytick,save_file):
    """Draw the biomass growth curve with different glucose uptake rates.
    Arguments:
    ----------
    * Biomass_list: The biomass growth curve with different glucose uptake rates.
    * save_file: The file name of the saved figure.
    """
    plt.figure(figsize=(15, 10), dpi=300)

    plt.plot(Biomass_list.index, Biomass_list[Biomass_list.columns[0]], color="black", linewidth=3.0, linestyle="--", label=Biomass_list.columns[0])
    plt.plot(Biomass_list.index, Biomass_list[Biomass_list.columns[1]], color="red", linewidth=3.0, linestyle="-", label=Biomass_list.columns[1])
    plt.plot(Biomass_list.index, Biomass_list[Biomass_list.columns[2]], color="cyan", linewidth=3.0, linestyle="-", label=Biomass_list.columns[2])
    plt.plot(Biomass_list.index, Biomass_list[Biomass_list.columns[3]], color="darkorange", linewidth=3.0, linestyle="-", label=Biomass_list.columns[3])

    font1 = {
    'weight' : 'normal',
    'size'   : 20,
    }

    plt.legend(loc="upper left",prop=font1)

    plt.xlim(xmin, xmax)
    plt.ylim(np.min(ytick), np.max(ytick))

    plt.tick_params(labelsize=23)
    plt.xticks(generate_list(xmin, xmax, step))
    plt.yticks(ytick)

    font2 = {
    'weight' : 'normal',
    'size'   : 25,
    }
    plt.xlabel("Glucose uptake rate (mmol/gDW/h)",font2)
    plt.ylabel("Growth rate ($\mathregular{h^-1}$)",font2)
    plt.savefig(save_file)
    plt.show()

def Draw_MDF_By_Growth_rate(MDF_list,xmin,xmax,step,ytick,save_file):
    """Draw the MDF curve with different growth rates.
    Arguments:
    ----------
    * MDF_list: The MDF curve with different growth rates.
    * save_file: The file name of the saved figure.
    """
    plt.figure(figsize=(15, 10), dpi=300)
    MDF_list=MDF_list.sort_index(ascending=True) 
    plt.plot(MDF_list.index, MDF_list[MDF_list.columns[0]], color="cyan", linewidth=3.0, linestyle="-", label=MDF_list.columns[0])
    plt.plot(MDF_list.index, MDF_list[MDF_list.columns[1]], color="darkorange", linewidth=3.0, linestyle="-", label=MDF_list.columns[1])
    font1 = {
    'weight' : 'normal',
    'size'   : 23,
    }

    font2 = {
    'weight' : 'normal',
    'size'   : 30,
    }

    plt.ylabel("MDF of pathways (kJ/mol)",font2)
    plt.xlabel("Growth rate ($\mathregular{h^-1}$)",font2)

    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    #ax.spines['left'].set_position(('data', 0.3))
    ax.spines['bottom'].set_position(('data', 0))
    plt.legend(loc="lower left",prop=font1)
    ax.spines['left'].set_position(('data', xmin))
    plt.xlim(xmin, xmax)
    plt.ylim(np.min(ytick), np.max(ytick))

    plt.tick_params(labelsize=23)

    plt.xticks(generate_list(xmin, xmax, step))
    plt.yticks(ytick)

    #plt.scatter([0.633], [2.6670879363966336], s=80, color="red")
    #plt.scatter([0.6756], [-0.48186643213771774], s=80, color="red")
    #plt.scatter([0.7068], [-9.486379882991386], s=80, color="red")
    #plt.scatter([0.852], [2.6670879363966336], s=80, color="red")
    #plt.scatter([0.855], [1.4290141211096987], s=80, color="red")
    #plt.scatter([0.867], [0.069495 15162540898], s=80, color="red")
    #plt.scatter([0.872], [-0.8364187795859692], s=80, color="red")
    #plt.scatter([0.876], [-9.486379882991372], s=80, color="red")

    plt.savefig(save_file)
    plt.show()
    
def Draw_MDF_By_Product_rate(MDF_dict, xmin, xlab, xmax, step, ytick, save_file):
    """Draw the MDF curves with different growth rates.
    Arguments:
    ----------
    * MDF_dict: A dictionary where keys are IDs and values are lists of points.
    * save_file: The file name of the saved figure.
    """
    plt.figure(figsize=(15, 10), dpi=300)

    font1 = {
        'weight': 'normal',
        'size': 23,
    }

    font2 = {
        'weight': 'normal',
        'size': 30,
    }

    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_position(('data', xmin))
    ax.spines['bottom'].set_position(('data', 0))

    for id, points in MDF_dict.items():
        x, y = zip(*points)
        plt.plot(x, y, label=id, linewidth=3.0)  # 使用不同的颜色

    plt.legend(loc="lower left", prop=font1)
    plt.ylabel("MDF of pathways (kJ/mol)", font2)
    plt.xlabel(xlab, font2)
    plt.xlim(xmin, xmax)
    plt.ylim(np.min(ytick), np.max(ytick))

    plt.tick_params(labelsize=23)
    plt.xticks(generate_list(xmin, xmax, step))
    plt.yticks(ytick)

    plt.savefig(save_file)
    plt.show()

def generate_list(min_value, max_value, step):
    result = []
    current_value = min_value
    while current_value <= max_value:
        result.append(current_value)
        current_value += step
    return result