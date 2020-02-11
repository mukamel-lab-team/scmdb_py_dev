from flask_assets import Bundle

app_css = Bundle('app_base.css', filters='cssmin', output='styles/app.css')

app_js = Bundle('app.js', filters='jsmin', output='scripts/app.js')

vendor_css = Bundle('vendor/semantic.min.css', output='styles/vendor.css')

vendor_js = Bundle(
    'vendor/jquery.min.js',
    'vendor/semantic.min.js',
    'vendor/tablesort.min.js',
    'vendor/zxcvbn.js',
    filters='jsmin',
    output='scripts/vendor.js')

browser_js = Bundle(
        'vendor/plotly.min.js',
        'vendor/datatables.min.js',
        'vendor/bootstrap-toggle.min.js',
        'vendor/select2.min.js',
        'vendor/bootstrap-slider.min.js',
        'customview.js',
        output='scripts/browser.js')
browser_css = Bundle(
        'vendor/datatables.min.css',
        'vendor/bootstrap-toggle.min.css',
        'vendor/select2.min.css',
        'vendor/bootstrap-slider.min.css',
        'browser.css',
        filters='cssmin',
        output='styles/browser.css')

tabular_rs1_js = Bundle(
        'vendor/datatables.min.js',
        'tabular_dataset_rs1.js',
        output='scripts/tabular_rs1.js')

tabular_rs2_js = Bundle(
        'vendor/datatables.min.js',
        'tabular_dataset_rs2.js',
        output='scripts/tabular_rs2.js')

tabular_ensemble_js = Bundle(
        'vendor/datatables.min.js',
        'tabular_ensemble.js',
        output='scripts/tabular_ensemble.js')

tabular_css = Bundle(
        'vendor/datatables.min.css',
        filters='cssmin',
        output='styles/tabular.css')

request_new_ensemble_js = Bundle(
        'vendor/datatables.min.js',
        'request_new_ensemble.js',
        output='scripts/request_new_ensemble.js')
