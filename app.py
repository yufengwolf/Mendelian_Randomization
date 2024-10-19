from flask import Flask, request, jsonify, render_template
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import json, os
import mr_database
import pandas as pd

# 启用 pandas 和 R 数据帧之间的转换
pandas2ri.activate()

# 设置 R 环境变量和库路径
robjects.r('''
source("~/.Rprofile")
''')

# 预先安装 R 包（仅在需要时运行一次）
def install_r_packages():
    robjects.r('''
    required_packages <- c("TwoSampleMR", "VariantAnnotation", "gwasvcf", "gwasglue", "ieugwasr", "dplyr", "ggplot2", "gridExtra", "forestploter")
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    }
    ''')
install_r_packages()  # 仅在需要时运行

# 延迟加载 R 包
def load_r_packages():
    robjects.r('''
    if (file.exists("~/.Renviron")) readRenviron("~/.Renviron")
    required_packages <- c("TwoSampleMR", "VariantAnnotation", "gwasvcf", "gwasglue", "ieugwasr", "dplyr", "ggplot2", "gridExtra", "forestploter")
    for (pkg in required_packages) {
      print(paste("Loading package:", pkg))
      library(pkg, character.only = TRUE)
    }
    ''')
    # 加载 get_opengwas_jwt 函数到全局环境
    robjects.r('get_opengwas_jwt <- ieugwasr::get_opengwas_jwt')

app = Flask(__name__)

exposure_dat = None
outcome_dat = None
harmonised_dat = None
mr_dat = None
het_dat = None
pleio_dat = None
leaveoneout_dat = None
leaveoneout_plot = None
singlesnp_dat = None
OR_dat = None
scatter_plot = None
forest_plot = None
funnel_plot = None
# 文件路径
excel_file = 'data.xlsx'

# 从文件中恢复数据
if os.path.exists(excel_file):
    with pd.ExcelFile(excel_file) as xls:
        if 'exposure_dat' in xls.sheet_names:
            exposure_dat = pd.read_excel(xls, sheet_name='exposure_dat')
        if 'outcome_dat' in xls.sheet_names:
            outcome_dat = pd.read_excel(xls, sheet_name='outcome_dat')
        if 'harmonised_dat' in xls.sheet_names:
            harmonised_dat = pd.read_excel(xls, sheet_name='harmonised_dat')
        if 'mr_dat' in xls.sheet_names:
            mr_dat = pd.read_excel(xls, sheet_name='mr_dat')
        if 'het_dat' in xls.sheet_names:
            het_dat = pd.read_excel(xls, sheet_name='het_dat')
        if 'pleio_dat' in xls.sheet_names:
            pleio_dat = pd.read_excel(xls, sheet_name='pleio_dat')
        if 'leaveoneout_dat' in xls.sheet_names:
            leaveoneout_dat = pd.read_excel(xls, sheet_name='leaveoneout_dat')
        if 'singlesnp_dat' in xls.sheet_names:
            singlesnp_dat = pd.read_excel(xls, sheet_name='singlesnp_dat')
        if 'OR_dat' in xls.sheet_names:
            OR_dat = pd.read_excel(xls, sheet_name='OR_dat')

def store_to_excel(data, sheet):
    if os.path.exists(excel_file):
        with pd.ExcelWriter(excel_file, mode='a', if_sheet_exists='replace') as writer:
            data.to_excel(writer, sheet_name=sheet, index=False)
    else:
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            data.to_excel(writer, sheet_name=sheet, index=False)

def extract_instruments():
    global exposure_dat  # 声明为全局变量
    outcome_id = request.form['outcome_id']
    p1 = float(request.form.get('p1', 5e-08))
    clump = request.form.get('clump', 'True') == 'True'
    p2 = float(request.form.get('p2', 5e-08))
    r2 = float(request.form.get('r2', 0.001))
    kb = int(request.form.get('kb', 10000))
    if not outcome_id:
        return None
    # 手动设置转换规则
    with localconverter(robjects.default_converter + pandas2ri.converter):
        extract_instruments = robjects.r['extract_instruments']
        get_opengwas_jwt = robjects.globalenv['get_opengwas_jwt']
        # print('get_opengwas_jwt: ', get_opengwas_jwt())
        exposure_dat = extract_instruments(outcomes=outcome_id, p1=p1, clump=clump, p2=p2, r2=r2, kb=kb, opengwas_jwt=get_opengwas_jwt(), force_server=False)
    store_to_excel(exposure_dat, 'exposure_dat')
    result = exposure_dat.to_json(orient='records')
    result = json.loads(result)
    return result

def extract_outcome():
    global exposure_dat, outcome_dat  # 声明为全局变量
    outcome_id = request.form['outcome_id']
    proxies = request.form.get('proxies', 'True') == 'True'
    rsq = float(request.form.get('rsq', 0.8))
    align_alleles = int(request.form.get('align_alleles', 1))
    palindromes = int(request.form.get('palindromes', 1))
    maf_threshold = float(request.form.get('maf_threshold', 0.3))
    splitsize = int(request.form.get('splitsize', 10000))
    proxy_splitsize = int(request.form.get('proxy_splitsize', 500))
    if not outcome_id:
        return None
    with localconverter(robjects.default_converter + pandas2ri.converter):
        extract_outcome_data = robjects.r['extract_outcome_data']
        get_opengwas_jwt = robjects.globalenv['get_opengwas_jwt']
        outcome_dat = extract_outcome_data(snps=exposure_dat['SNP'], outcomes=outcome_id, proxies=proxies, rsq=rsq, align_alleles=align_alleles, palindromes=palindromes, maf_threshold=maf_threshold, opengwas_jwt=get_opengwas_jwt(), splitsize=splitsize, proxy_splitsize=proxy_splitsize)
        store_to_excel(outcome_dat, 'outcome_dat')
    result = outcome_dat.to_json(orient='records')
    result = json.loads(result)
    return result

def extract_harmonise():
    global exposure_dat, outcome_dat, harmonised_dat  # 声明为全局变量
    action = int(request.form.get('action', 2))
    with localconverter(robjects.default_converter + pandas2ri.converter):
        harmonise_data = robjects.r['harmonise_data']
        harmonised_dat = harmonise_data(exposure_dat=exposure_dat, outcome_dat=outcome_dat, action=action)
        store_to_excel(harmonised_dat, 'harmonised_dat')
    result = harmonised_dat.to_json(orient='records')
    result = json.loads(result)
    return result

def mr_analysis():
    global harmonised_dat, mr_dat
    method_list = request.form['method_list']
    if not method_list:
        return None
    with localconverter(robjects.default_converter + pandas2ri.converter):
        mr = robjects.r['mr']
        c = robjects.r['c']
        mr_dat = mr(harmonised_dat, method_list=c(method_list))
        store_to_excel(mr_dat, 'mr_dat')
    result = mr_dat.to_json(orient='records')
    result = json.loads(result)
    return result

def mr_heterogeneity_analysis():
    global harmonised_dat, het_dat
    with localconverter(robjects.default_converter + pandas2ri.converter):
        mr_heterogeneity = robjects.r['mr_heterogeneity']
        het_dat = mr_heterogeneity(harmonised_dat)
        store_to_excel(het_dat, 'het_dat')
    result = het_dat.to_json(orient='records')
    result = json.loads(result)
    return result

def mr_pleiotropy_test_analysis():
    global harmonised_dat, pleio_dat
    # 手动设置转换规则
    with localconverter(robjects.default_converter + pandas2ri.converter):
        pleiotropy_test = robjects.r['mr_pleiotropy_test']
        pleio_dat = pleiotropy_test(harmonised_dat)
        store_to_excel(pleio_dat, 'pleio_dat')
    result = pleio_dat.to_json(orient='records')
    result = json.loads(result)
    return result

def mr_leaveoneout_analysis():
    global harmonised_dat, leaveoneout_dat, leaveoneout_plot
   # 手动设置转换规则
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_harmonised_dat = pandas2ri.py2rpy(harmonised_dat)
        robjects.globalenv['harmonised_dat'] = r_harmonised_dat
        plot_dir = os.path.abspath('static')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        leaveoneout_plot = 'leaveoneout_plot.png'
        leaveoneout_plot_path = os.path.join(plot_dir, leaveoneout_plot)
        robjects.globalenv['leaveoneout_plot_path'] = leaveoneout_plot_path

        r_code = '''
            leaveoneout_dat <- mr_leaveoneout(harmonised_dat)
            leaveoneout_plot <- mr_leaveoneout_plot(leaveoneout_dat)
            ggsave(leaveoneout_plot_path, leaveoneout_plot[[1]])
        '''
        # print("R code to be executed:", r_code)  # Debugging information
        robjects.r(r_code)
        leaveoneout_dat = robjects.globalenv['leaveoneout_dat']
        store_to_excel(leaveoneout_dat, 'leaveoneout_dat')

    result = leaveoneout_dat.to_json(orient='records')
    result = json.loads(result)
    return result, leaveoneout_plot

def mr_singlesnp_analysis():
    global harmonised_dat, singlesnp_dat
    # 手动设置转换规则
    with localconverter(robjects.default_converter + pandas2ri.converter):
        mr_singlesnp = robjects.r['mr_singlesnp']
        singlesnp_dat = mr_singlesnp(harmonised_dat)
        store_to_excel(singlesnp_dat, 'singlesnp_dat')
    result = singlesnp_dat.to_json(orient='records')
    result = json.loads(result)
    return result

def generate_odds_ratios():
    global mr_dat, OR_dat
    # 手动设置转换规则
    with localconverter(robjects.default_converter + pandas2ri.converter):
        generate_odds_ratios = robjects.r['generate_odds_ratios']
        OR_dat = generate_odds_ratios(mr_dat)
        store_to_excel(OR_dat, 'OR_dat')
    result = OR_dat.to_json(orient='records')
    result = json.loads(result)
    return result

def mr_scatter_plot():
    global mr_dat, harmonised_dat, scatter_plot
    # 手动设置转换规则
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_harmonised_dat = pandas2ri.py2rpy(harmonised_dat)
        r_mr_dat = pandas2ri.py2rpy(mr_dat)
        robjects.globalenv['harmonised_dat'] = r_harmonised_dat
        robjects.globalenv['mr_dat'] = r_mr_dat
        plot_dir = os.path.abspath('static')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        scatter_plot = 'scatter_plot.png'
        scatter_plot_path = os.path.join(plot_dir, scatter_plot)
        robjects.globalenv['scatter_plot_path'] = scatter_plot_path
        r_code = '''
            scatter_plot <- mr_scatter_plot(mr_dat, harmonised_dat)
            ggsave(scatter_plot_path, scatter_plot[[1]])
        '''
        robjects.r(r_code)
    return None, scatter_plot

def mr_forest_plot():
    global singlesnp_dat, forest_plot
    # 手动设置转换规则
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_singlesnp_dat = pandas2ri.py2rpy(singlesnp_dat)
        robjects.globalenv['singlesnp_dat'] = r_singlesnp_dat
        plot_dir = os.path.abspath('static')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        forest_plot = 'forest_plot.png'
        forest_plot_path = os.path.join(plot_dir, forest_plot)
        robjects.globalenv['forest_plot_path'] = forest_plot_path
        r_code = '''
            forest_plot <- mr_forest_plot(singlesnp_dat)
            ggsave(forest_plot_path, forest_plot[[1]])
        '''
        robjects.r(r_code)
    return None, forest_plot

def mr_funnel_plot():
    global singlesnp_dat, funnel_plot
    # 手动设置转换规则
    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_singlesnp_dat = pandas2ri.py2rpy(singlesnp_dat)
        robjects.globalenv['singlesnp_dat'] = r_singlesnp_dat
        plot_dir = os.path.abspath('static')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        funnel_plot = 'funnel_plot.png'
        funnel_plot_path = os.path.join(plot_dir, funnel_plot)
        robjects.globalenv['funnel_plot_path'] = funnel_plot_path
        r_code = '''
            funnel_plot <- mr_funnel_plot(singlesnp_dat)
            ggsave(funnel_plot_path, funnel_plot[[1]])
        '''
        robjects.r(r_code)
    return None, funnel_plot

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/show/exposure_data')
def show_exposure_data():
    global exposure_dat
    result = exposure_dat.to_json(orient='records')
    return render_template('result.html', result=json.loads(result), plot_path='')

@app.route('/show/outcome_data')
def show_outcome_data():
    global outcome_dat
    result = outcome_dat.to_json(orient='records')
    return render_template('result.html', result=json.loads(result), plot_path='')

@app.route('/show/harmonised_data')
def show_harmonised_data():
    global harmonised_dat
    result = harmonised_dat.to_json(orient='records')
    return render_template('result.html', result=json.loads(result), plot_path='')

@app.route('/show/mr_data')
def show_mr_data():
    global mr_dat
    result = mr_dat.to_json(orient='records')
    return render_template('result.html', result=json.loads(result), plot_path='')

@app.route('/show/OR_data')
def show_OR_data():
    global OR_dat
    result = OR_dat.to_json(orient='records')
    return render_template('result.html', result=json.loads(result), plot_path='')

@app.route('/show/heterogeneity_data')
def show_het_data():
    global het_dat
    result = het_dat.to_json(orient='records')
    return render_template('result.html', result=json.loads(result), plot_path='')

@app.route('/show/pleiotropy_test_data')
def show_pleio_data():
    global pleio_dat
    result = pleio_dat.to_json(orient='records')
    return render_template('result.html', result=json.loads(result), plot_path='')

@app.route('/show/leaveoneout_data')
def show_leaveoneout_data():
    global leaveoneout_dat, leaveoneout_plot
    result = leaveoneout_dat.to_json(orient='records')
    return render_template('result.html', result=json.loads(result), plot_path=leaveoneout_plot)

@app.route('/show/singlesnp_data')
def show_singlesnp_data():
    global singlesnp_dat
    result = singlesnp_dat.to_json(orient='records')
    return render_template('result.html', result=json.loads(result), plot_path='')

@app.route('/show/scatter_plot')
def show_scatter_plot():
    global scatter_plot
    return render_template('result.html', result=None, plot_path=scatter_plot)

@app.route('/show/forest_plot')
def show_forest_plot():
    global forest_plot
    return render_template('result.html', result=None, plot_path=forest_plot)

@app.route('/show/funnel_plot')
def show_funnel_plot():
    global funnel_plot
    return render_template('result.html', result=None, plot_path=funnel_plot)

@app.route('/call_function', methods=['POST'])
def call_function():
    global exposure_dat, outcome_dat, harmonised_dat, mr_dat  # 声明为全局变量
    load_r_packages()  # 延迟加载 R 包
    function = request.form['function']
    result = None
    plot_path = ''
    if function == 'extract_instruments':
        result = extract_instruments()
    elif function == 'extract_outcome':
        result = extract_outcome()
    elif function == 'extract_harmonise':
        result = extract_harmonise()
    elif function == 'mr_analysis':
        result = mr_analysis()
    elif function == 'generate_odds_ratios':
        result = generate_odds_ratios()
    elif function == 'mr_heterogeneity_analysis':
        result = mr_heterogeneity_analysis()
    elif function == 'mr_pleiotropy_test_analysis':
        result = mr_pleiotropy_test_analysis()
    elif function == 'mr_leaveoneout_analysis':
        result, plot_path = mr_leaveoneout_analysis()
    elif function == 'mr_singlesnp_analysis':
        result = mr_singlesnp_analysis()
    elif function == 'mr_scatter_plot':
        result, plot_path = mr_scatter_plot()
    elif function == 'mr_forest_plot':
        result, plot_path = mr_forest_plot()
    elif function == 'mr_funnel_plot':
        result, plot_path = mr_funnel_plot()
    else:
        return jsonify({'error': 'Invalid function'}), 400
    if result == None:
        jsonify({'error': 'Missing parameter'}), 400
    return render_template('result.html', result=result, plot_path=plot_path)

if __name__ == '__main__':
    app.run(debug=True)