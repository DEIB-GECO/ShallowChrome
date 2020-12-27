# python auxiliary routines

from sklearn import metrics as met
import pickle
import numpy as np
from scipy.stats import norm, zscore
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.cm import get_cmap

def load_cell_paths(cell_list_path, data_base_path, suffix='', best_formats=None, cells=None):
    
    if cells is None:
        cells = np.loadtxt(cell_list_path, dtype=str, delimiter=',')
    try:
        cells = cells.reshape((len(cells),))
    except TypeError:
        cells = cells.reshape((1,))
    
    if best_formats is not None:
        folder_list = [data_base_path+best_formats[cell]+'/'+cell+suffix for cell in cells]
    else:
        folder_list = [data_base_path+cell+suffix for cell in cells]

    return folder_list, cells


def binarize_target(T, threshold=None, percentile=50):

    if threshold is None:
        threshold = np.percentile(T, percentile)

    T_bin = np.zeros(T.shape)
    T_bin[T>=threshold] = 1

    return T_bin


def load_split(base_split, suffix_split, split_idx):

    folder = base_split+suffix_split+str(split_idx)+'/'

    _train = np.loadtxt(folder+'_train.csv', delimiter=',', dtype=np.int32)
    _val = np.loadtxt(folder+'_val.csv', delimiter=',', dtype=np.int32)
    _test = np.loadtxt(folder+'_test.csv', delimiter=',', dtype=np.int32)

    return (_train, _val, _test)


def random_splits(n, iterations=20, cache_at=None, random_state=None):

    random_obj = np.random.RandomState(random_state)
    step = n/3
    train = [i for i in range(step+1)]
    val = [i for i in range(step+1,2*(step+1))]
    test = [i for i in range(2*(step+1),n)]

    if cache_at is not None:
        try:
            os.mkdir(cache_at)
        except OSError:
            pass

    splits = []
    for i in range(iterations):

        rand = random_obj.permutation(n)
        splits.append((rand[train], rand[val], rand[test]))

        if cache_at is not None:
            iteration_path = cache_at+'/iteration_'+str(i+1)+'/'
            try:
                os.mkdir(iteration_path)
            except OSError as e:
                pass
            np.savetxt(iteration_path+'_train.csv', rand[train], delimiter=',', fmt='%d')
            np.savetxt(iteration_path+'_val.csv', rand[val], delimiter=',', fmt='%d')
            np.savetxt(iteration_path+'_test.csv', rand[test], delimiter=',', fmt='%d')

    return splits


def balance_by_subsampling(X, T, random_state=None):

    n_pos = (T==1).sum()
    n_neg = len(T) - n_pos
    pos_samples = np.where(T==1)
    neg_samples = np.where(T==0)
    pos_rand = np.random.RandomState(random_state).permutation(n_pos)
    neg_rand = np.random.RandomState(random_state).permutation(n_neg)

    if n_pos > n_neg:
        target = n_neg
    elif n_neg > n_pos:
        target = n_pos
    else:
        return (X, T)

    X_pos_sub = X[pos_samples][pos_rand][:target]
    X_neg_sub = X[neg_samples][neg_rand][:target]
    T_pos_sub = T[pos_samples][pos_rand][:target]
    T_neg_sub = T[neg_samples][neg_rand][:target]

    X_sub = np.concatenate((X_pos_sub, X_neg_sub), axis=0)
    T_sub = np.concatenate((T_pos_sub, T_neg_sub), axis=0)

    return (X_sub, T_sub)



def fit_and_score(model, X, T, split, cache_model_at=None):

    _train, _val, _test = split

    model.fit(X[_train], T[_train])
    if cache_model_at is not None:
        with open(cache_model_at, 'wb') as model_file:
            pickle.dump(model, model_file)

    val_score = met.roc_auc_score(T[_val], model.decision_function(X[_val]))
    test_score = met.roc_auc_score(T[_test], model.decision_function(X[_test]))

    return (val_score, test_score)


def z_test(X, y, model, names, alpha=None):
    n_samples, n_features = X.shape
    deg = n_samples-n_features
    betas = np.concatenate([model.intercept_, model.coef_.reshape(-1)])
    
    # Compute the prediction
    pred = model.predict_proba(X) # [N, 2]
    y = y.reshape(-1)    
    X = np.concatenate([np.ones([X.shape[0], 1]), X], axis=-1)
    n_samples, n_features = X.shape
    
    V = np.diagflat(np.product(pred, axis=1))
    covLogit = np.linalg.inv(np.dot(np.dot(X.T, V), X))
    se_b = np.sqrt(np.diag(covLogit))
    
    z_stat_b = (betas - 0) / se_b

    # Compute the p-value (two-sided test)
    p_values = np.array([2 * norm.sf(np.abs(z_stat)) for z_stat in z_stat_b])
    
    df = pd.DataFrame()
    df["Name"] = names
    df["Coefficients"] = betas
    df["Standard Errors"] = np.round(se_b, decimals=4)
    df["Z-stat"] = np.round(z_stat_b, decimals=1)
    df["p-value"] = p_values
    if alpha:
        rejectH0 = p_values < alpha
        df["reject H0"] = rejectH0    
    
    return df


def extract_patterns(model, X, collapse=True):

    w = model.coef_
    b = model.intercept_[0]
    
    if collapse:
        pattern = np.asarray([b]+np.median(X*w, axis=0).tolist())
        up = np.asarray([b]+np.percentile(X*w, 75, axis=0).tolist())
        down = np.asarray([b]+np.percentile(X*w, 25, axis=0).tolist())
        return (pattern, up, down)
    
    else:
        biases = np.asarray([b]*len(X)).reshape((len(X), 1))
        psis = X*w
        patterns = np.hstack((biases, psis))
        return patterns, b, psis
    
    
def draw_specific_pattern(
        model,
        pattern,
        name_array,
        X_train,
        y_train,
        show_intercept=True,
        absolute=True,
        norm=True,
        y_bounds=None,
        cache_at=None,
        label=None,
        legend=True,
        show=False,
        alpha=0.0001):
    
    # inspect significance
    df = z_test(X_train, y_train, model, ['(bias)'] + name_array.tolist(), alpha=alpha)
    keep = np.asarray(df['reject H0'])[1:]
    
    # get roles
    w = model.coef_
    repressors = np.where(np.logical_and(w<0, keep))
    activators = np.where(np.logical_and(w>0, keep))
    discard = np.where(~keep)
    
    # intercept and absolute?
    if show_intercept:
        pattern = np.concatenate((model.intercept_[np.newaxis, :], pattern), axis=1)
    if absolute:
        pattern = abs(pattern)
    if norm:
        pattern = pattern / np.sum(abs(pattern))

    # x positions
    x = np.arange(1, len(pattern[0])+1, 1)

    # y bounds
    if y_bounds is None:
        y_bounds = [min(pattern[0])-0.05, max(pattern[0])+0.05]

    # plot
    fig = plt.figure(dpi=300, figsize=(6, 4))
    plt.ylim(y_bounds)
    if show_intercept:
        barlist = plt.bar(x, pattern[0], width=0.8)
        for r in repressors[-1]:
            barlist[r+1].set_color('indianred')
        for a in activators[-1]:
            barlist[a+1].set_color('cadetblue')
        for d in discard[-1]:
            barlist[d+1].set_color('lightgray')
    else:
        barlist = plt.bar(x, pattern[0], width=0.8)
        for r in repressors[-1]:
            barlist[r].set_color('indianred')
        for a in activators[-1]:
            barlist[a].set_color('cadetblue')
        for d in discard[-1]:
            barlist[d+1].set_color('lightgray')
    
    if show_intercept:
        color = 'indianred' if model.intercept_ <= 0 else 'cadetblue'
        barlist[0].set_color(color)
    barlist = plt.bar(x, np.zeros(x.shape), width=0.8, color='cadetblue', label='activator')
    barlist = plt.bar(x, np.zeros(x.shape), width=0.8, color='indianred', label='repressor')
    
    if legend:
        plt.legend(fontsize=10)
                
    # manage tick labels and title
    if show_intercept:
        labels = np.hstack(('(bias)', name_array))
    else:
        labels = name_array
    plt.xticks(x, labels.tolist(), rotation=90)
    y_ticks = [ min([y_bounds[-1],pattern[0][0]]), pattern[0][1:].min(), min([y_bounds[-1],pattern[0][1:].max()]) ]
    y_labels = ['{0:.1f} %'.format(t*100) for t in y_ticks]
    if show_intercept:
        y_labels[0] = ''+'{0:.1f} %'.format(pattern[0][0]*100)+''
    plt.yticks(y_ticks, y_labels, fontsize=11)
    tail = ''
    if label is not None:
        tail = str(label)
    plt.title(''+tail)

    # cache
    if cache_at is not None:
        plt.savefig(cache_at, dpi=300, format='pdf', bbox_inches='tight', pad_inches=.0)
    if show:
        plt.show()
    plt.close()

    return

    
def extract_emissions(emission_path, delimiter='\t'):
    
    hm_index = dict()
    hs_index = dict()
    emissions = list()
    
    with open(emission_path, 'r') as emi_file:
        
        header = next(emi_file).strip().split(delimiter)[1:]
        for h, head in enumerate(header):
            hm_index[h] = head
        
        for l, line in enumerate(emi_file):
            fields = line.strip().split(delimiter)
            hs_index[l] = fields[0]
            float_fields = [float(field) for field in fields[1:]]
            emissions.append(np.asarray(float_fields))
            
        emissions = np.asarray(emissions)
    
    return emissions, hm_index, hs_index


def extract_states_colors(state_path, delimiter='\t', state_name=1, state_description=2, color_code=4):
    
    index = dict()
    color_dict = dict()
    description_dict = dict()
    with open(state_path, 'r') as s_file:
        
        header = next(s_file).strip().split(delimiter)
        
        for l, line in enumerate(s_file):
            if l==0:
                continue
            fields = line.strip().split(delimiter)
            name = fields[state_name]
            descr = fields[state_description]
            color = [float(channel)/255.0 for channel in fields[color_code].split(',')]
            index[l] = name
            description_dict[name] = descr
            color_dict[name] = color
    
    return index, description_dict, color_dict


def compute_mean_std_logits_per_state(logits, bmus, num_states):
    
    mean_logits_by_state = dict()
    std_logits_by_state = dict()
    for state in range(num_states):
        indexes = np.where(bmus==state)
        matching_logits = logits[indexes]
        if len(matching_logits)==0:
            mean_logits_by_state[state+1] = -np.inf
            std_logits_by_state[state+1] = -np.inf
        else:
            mean_logits_by_state[state+1] = np.mean(matching_logits)
            std_logits_by_state[state+1] = np.std(matching_logits)
    mean_logits = np.asarray([mean_logits_by_state[state] for state in range(1, num_states+1)])
    std_logits = np.asarray([std_logits_by_state[state] for state in range(1, num_states+1)])
    
    return mean_logits, std_logits


def find_valley(T_tr, freqs, edges):
    
    # define proper interval
    end = np.percentile(T_tr, 50, interpolation='linear')
    slicer = np.where(edges<end)
    
    # find valley as minimum point in there
    minimum = np.argmin(freqs[slicer])
    valley_tr = edges[minimum]
    
    return valley_tr


def parse_scores(score_path, std=False):

    scores = {}
    with open(score_path, 'r') as score_file:
        for line in score_file:
            fields = line.strip().split(': ')
            if std is False:
                scores[fields[0]] = float(fields[1])
            else:
                values = fields[1].split(' +/- ')
                scores[fields[0]] = (float(values[0]), float(values[1]))

    return _extract_score_vector(scores, std=std)


def _extract_score_vector(score_dict, std=False):

    cells = list(score_dict.keys())
    cells.sort()
    if std is False:
        score_vector = np.ndarray((len(cells), 1))
        for c, cell in enumerate(cells):
            score_vector[c] = score_dict[cell]
    else:
        score_vector = np.ndarray((len(cells), 2))
        for c, cell in enumerate(cells):
            score_vector[c,0], score_vector[c,1] = score_dict[cell][0], score_dict[cell][1]

    return score_vector
