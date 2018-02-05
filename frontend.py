"""Defines the web routes of the website.

For actual content generation see the content.py module.
"""
from os import walk
import os.path

import dominate
from dominate.tags import img
from flask import Blueprint, render_template, jsonify, request, redirect, current_app, flash, abort, url_for
from flask_login import (current_user, login_required, login_user,
                         logout_user)
from flask_mail import Mail, Message
from flask_nav.elements import Navbar, Link, View, Text
from flask_rq import get_queue

from . import nav, cache, db, mail
from .content import get_cluster_plot, search_gene_names, \
    get_methylation_scatter, get_mch_box, get_mch_box_two_dataset, \
    find_orthologs, FailToGraphException, get_corr_genes, \
    gene_id_to_name, randomize_cluster_colors, get_mch_heatmap, \
    get_mch_heatmap_two_dataset, all_gene_modules, get_genes_of_module, \
    convert_geneID_mmu_hsa
from .decorators import admin_required
from .email import send_email
from .forms import LoginForm, ChangeUserEmailForm, ChangeAccountTypeForm, InviteUserForm, CreatePasswordForm, NewUserForm, RequestResetPasswordForm, ResetPasswordForm, ChangePasswordForm
from .user import User, Role


frontend = Blueprint('frontend', __name__) # Flask "bootstrap"


# HOTFIX: '/' character needed to prevent concatenation of url
@frontend.before_request
def process_navbar():
    # get images here
    lockimage = img(src='static/img/lock.png', height='20', width='20')
    unlockimage = img(src='static/img/unlock.png', height='20', width='20')
    separator = img(src='static/img/separate.png', height='25', width='10')

    # Find all the samples in the data directory
    dir_list = []
    if current_user.is_authenticated:
        current_app.config['DATA_DIR'] = current_app.config['ALL_DATA_DIR']
    else:
        current_app.config['DATA_DIR'] = current_app.config['PUBLISHED_DATA_DIR']
    dir_list = next(walk(current_app.config['DATA_DIR'] + "/ensembles"))[1]
    published_dir_list = next(walk(current_app.config['PUBLISHED_DATA_DIR'] + 
            "/ensembles"))[1]

    dir_list_links = []

    first = True

    for x in dir_list:
        if not first:
            dir_list_links.append(Text(separator))
        dir_list_links.append(Link(x, "/" + x))
        if current_user.is_authenticated:
            # x is public: add unlockimage
            if x in published_dir_list:
                dir_list_links.append(Text(unlockimage))
            # x is private: add lockimage
            else:
                dir_list_links.append(Text(lockimage))
        first = False

    nav.register_element('frontend_top', Navbar('',*dir_list_links))

# Visitor routes
@frontend.route('/')
def index():
    # Index is not needed since this site is embedded as a frame
    # We use JS redirect b/c reverse proxy will be confused about subdirectories
    html = \
    """To be redirected manually, click <a href="./hsa">here</a>.
    <script>
        window.location = "./human_hv1_published"; 
        window.location.replace("./human_hv1_published");
    </script>
    """
    return html
    
    # TODO: May need to switch to code below 
    #    depending on how to deal w/ published data

    # return \
    # """
    # To be redirected manually, click <a href="./human_combined">here</a>.' + \
    # <script>
    #      window.location = "./human_combined"; 
    #      window.location.replace("./human_combined");
    # </script>
    # """

@frontend.route('/<dataset>')
def dataset(dataset):
    return render_template('datasetview.html', dataset=dataset)

@frontend.route('/standalone/<dataset>/<gene>')
def standalone(dataset, gene):  # View gene body mCH plots alone
    return render_template('mch_standalone.html', dataset=dataset, gene=gene)


@frontend.route('/compare/<mmu_gid>/<hsa_gid>')
def compare(mmu_gid, hsa_gid):
    return render_template('compareview.html', mmu_gid=mmu_gid, hsa_gid=hsa_gid)


@frontend.route('/box_combined/<mmu_gid>/<hsa_gid>')
def box_combined(mmu_gid, hsa_gid):
    return render_template(
        'combined_box_standalone.html', mmu_gid=mmu_gid, hsa_gid=hsa_gid)


@frontend.route('/tabular/ensemble')
def ensemble_tabular_screen():
    return render_template('tabular_ensemble.html')


@frontend.route('/tabular/dataset')
def data_set_tabular_screen():
    return render_template('tabular_data_set.html')


@frontend.route('/navbar')
def nav_bar_screen():
    return render_template('navbar_only.html')


# API routes
@frontend.route('/plot/cluster/<dataset>/<grouping>')
@cache.memoize(timeout=3600)
def plot_cluster(dataset, grouping):
    try:
        return jsonify(get_cluster_plot(dataset, grouping))
    except FailToGraphException:
        return 'Failed to produce cluster plot. Contact maintainer.'


@frontend.route('/plot/scatter/<dataset>/<methylationType>/<level>/<ptile_start>/<ptile_end>')
def plot_methylation_scatter(dataset, methylationType, level, ptile_start, ptile_end):
    genes = request.args.get('q', 'MustHaveAQueryString')
    try:
        return get_methylation_scatter(dataset,
                                       methylationType,
                                       genes, level,
                                       float(ptile_start), float(ptile_end))
    except (FailToGraphException, ValueError) as e:
        print(e)
        return 'Failed to produce mCH levels scatter plot. Contact maintainer.'


@frontend.route('/plot/box/<dataset>/<methylationType>/<gene>/<level>/<outliers_toggle>')
@cache.memoize(timeout=3600)
def plot_mch_box(dataset, methylationType, gene, level, outliers_toggle):
    if outliers_toggle == 'outliers':
        outliers = True
    else:
        outliers = False

    try:
        return get_mch_box(dataset, methylationType, gene, level, outliers)
    except (FailToGraphException, ValueError) as e:
        print(e)
        return 'Failed to produce mCH levels box plot. Contact maintainer.'


@frontend.route('/plot/box_combined/<methylationType>/<gene_mmu>/<gene_hsa>/<level>/<outliers_toggle>')
def plot_mch_box_two_dataset(methylationType, gene_mmu, gene_hsa, level, outliers_toggle):
    if outliers_toggle == 'outliers':
        outliers = True
    else:
        outliers = False
    try:
        return get_mch_box_two_dataset(methylationType, gene_mmu, gene_hsa, level, outliers)
    except (FailToGraphException, ValueError) as e:
        print(e)
        return 'Failed to produce mCH levels box plot. Contact maintainer.'


@frontend.route('/plot/heat/<dataset>/<methylationType>/<level>/<ptile_start>/<ptile_end>')
def plot_mch_heatmap(dataset, methylationType, level, ptile_start, ptile_end):
    query = request.args.get('q', 'MustHaveAQueryString')
    if request.args.get('normalize', 'MustSpecifyNormalization') == 'true':
        normalize_row = True
    else:
        normalize_row = False
    try:
        return get_mch_heatmap(dataset, methylationType, level, ptile_start, ptile_end, normalize_row, query)
    except (FailToGraphException, ValueError) as e:
        print(e)
        return 'Failed to produce mCH levels heatmap plot. Contact maintainer.'

@frontend.route('/plot/heat_two_dataset/<dataset>/<methylationType>/<level>/<ptile_start>/<ptile_end>')
def plot_mch_heatmap_two_dataset(dataset, methylationType, level, ptile_start, ptile_end):
    query = request.args.get('q', 'MustHaveAQueryString')
    if request.args.get('normalize', 'MustSpecifyNormalization') == 'true':
        normalize_row = True
    else:
        normalize_row = False
    try:
        return get_mch_heatmap_two_dataset(dataset, methylationType, level, ptile_start, ptile_end, normalize_row, query)
    except (FailToGraphException, ValueError) as e:
        print(e)
        return 'Failed to produce orthologous mCH levels heatmap plot. Contact maintainer.'

@frontend.route('/gene/names/<dataset>')
def search_gene_by_name(dataset):
    query = request.args.get('q', 'MustHaveAQueryString')
    if query == 'none' or query == '':
        return jsonify([])
    else:
        return jsonify(search_gene_names(dataset, query))


@frontend.route('/gene/id/<dataset>')
def search_gene_by_id(dataset):
    query = request.args.get('q', 'MustHaveAQueryString')
    if query == 'none' or query == '':
        return jsonify({})
    else:
        return jsonify(gene_id_to_name(dataset, convert_geneID_mmu_hsa(dataset, query)))


@frontend.route('/gene/modules/<dataset>')
def gene_modules(dataset):
    query = request.args.get('q')
    if query == None or query == '':
        return jsonify(all_gene_modules())
    else:
        return jsonify(get_genes_of_module(dataset, query))


@frontend.route('/gene/orthologs/<dataset>/<geneID>')
def orthologs(dataset, geneID):
    geneID = geneID.split('.')[0]
    geneID = convert_geneID_mmu_hsa(dataset,geneID) 
    if dataset == 'mmu' or dataset == 'mouse_published':
        return jsonify(find_orthologs(mmu_gid=geneID))
    else:
        return jsonify(find_orthologs(hsa_gid=geneID))


@cache.memoize(timeout=3600)
@frontend.route('/gene/corr/<dataset>/<geneID>')
def correlated_genes(dataset, geneID):
    return jsonify(get_corr_genes(dataset, geneID))


@frontend.route('/plot/randomize_colors')
def randomize_colors():
    num_colors = request.args.get('n', type=int)
    return jsonify(randomize_cluster_colors(num_colors))


@frontend.route('/plot/delete_cache/<dataset>/<grouping>')
def delete_cluster_cache(dataset, grouping):
    cache.delete_memoized(plot_cluster, dataset, grouping)
    return (dataset + " cluster cache cleared") 


# User related routes
@frontend.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        if user is not None and user.password_hash is not None and \
                user.verify_password(form.password.data):
            login_user(user, form.remember_me.data)
            return redirect('/')
        else:
            flash('Invalid email or password.', 'form-error')
    return render_template('account/login.html', form=form)

@frontend.route('/logout')
@login_required
def logout():
    logout_user()
    flash('You have been logged out.', 'info')
    return redirect(url_for('frontend.index'))
  
@frontend.route('/admin')
@login_required
@admin_required
def admin():
    """Admin dashboard page."""
    return render_template('admin/index.html')

@frontend.route('/users')
@login_required
@admin_required
def registered_users():
    """View all registered users."""
    users = User.query.all()
    roles = Role.query.all()
    return render_template(
        'admin/registered_users.html', users=users, roles=roles)

@frontend.route('/user/<int:user_id>')
@frontend.route('/user/<int:user_id>/info')
@login_required
@admin_required
def user_info(user_id):
    """View a user's profile."""
    user = User.query.filter_by(id=user_id).first()
    if user is None:
        abort(404)
    return render_template('admin/manage_user.html', user=user)


@frontend.route('/user/<int:user_id>/change-email', methods=['GET', 'POST'])
@login_required
@admin_required
def change_user_email(user_id):
    """Change a user's email."""
    user = User.query.filter_by(id=user_id).first()
    if user is None:
        abort(404)
    form = ChangeUserEmailForm()
    if form.validate_on_submit():
        user.email = form.email.data
        db.session.add(user)
        db.session.commit()
        flash('Email for user {} successfully changed to {}.'
              .format(user.full_name(), user.email), 'form-success')
    return render_template('admin/manage_user.html', user=user, form=form)

@frontend.route('/user/<int:user_id>/delete')
@login_required
@admin_required
def delete_user_request(user_id):
    """Request deletion of a user's account."""
    user = User.query.filter_by(id=user_id).first()
    if user is None:
        abort(404)
    return render_template('admin/manage_user.html', user=user)

@frontend.route('/user/<int:user_id>/_delete')
@login_required
@admin_required
def delete_user(user_id):
    """Delete a user's account."""
    if current_user.id == user_id:
        flash('You cannot delete your own account. Please ask another '
              'administrator to do this.', 'error')
    else:
        user = User.query.filter_by(id=user_id).first()
        db.session.delete(user)
        db.session.commit()
        flash('Successfully deleted user %s.' % user.full_name(), 'success')
    return redirect(url_for('frontend.registered_users'))


@frontend.route('/user/<int:user_id>/change-account-type', methods=['GET', 'POST'])
@login_required
@admin_required
def change_account_type(user_id):
    """Change a user's account type."""
    if current_user.id == user_id:
        flash('You cannot change the type of your own account. Please ask '
              'another administrator to do this.', 'error')
        return redirect(url_for('frontend.user_info', user_id=user_id))

    user = User.query.get(user_id)
    if user is None:
        abort(404)
    form = ChangeAccountTypeForm()
    if form.validate_on_submit():
        user.role = form.role.data
        db.session.add(user)
        db.session.commit()
        flash('Role for user {} successfully changed to {}.'
              .format(user.full_name(), user.role.name), 'form-success')
    return render_template('admin/manage_user.html', user=user, form=form)


@frontend.route('/invite-user', methods=['GET', 'POST'])
@login_required
@admin_required
def invite_user():
    """Invites a new user to create an account and set their own password."""
    form = InviteUserForm()
    if form.validate_on_submit():
        user = User(
            role=form.role.data,
            first_name=form.first_name.data,
            last_name=form.last_name.data,
            email=form.email.data)
        db.session.add(user)
        db.session.commit()
        token = user.generate_confirmation_token()
        invite_link = url_for(
            'frontend.join_from_invite',
            user_id=user.id,
            token=token,
            _external=True)
        send_email(recipient=user.email,subject='You Are Invited To Join',template='email/invite',user=user,invite_link=invite_link)
        flash('User {} successfully invited'.format(user.full_name()),
              'form-success')
    return render_template('admin/new_user.html', form=form)


@frontend.route(
    '/join-from-invite/<int:user_id>/<token>', methods=['GET', 'POST'])
def join_from_invite(user_id, token):
    """
    Confirm new user's account with provided token and prompt them to set
    a password.
    """
    if current_user is not None and current_user.is_authenticated:
        flash('You are already logged in.', 'error')
        return redirect(url_for('frontend.index'))

    new_user = User.query.get(user_id)
    if new_user is None:
        return redirect(404)

    if new_user.password_hash is not None:
        flash('You have already joined.', 'error')
        return redirect(url_for('frontend.index'))

    if new_user.confirm_account(token):
        form = CreatePasswordForm()
        if form.validate_on_submit():
            new_user.password = form.password.data
            db.session.add(new_user)
            db.session.commit()
            flash('Your password has been set. After you log in, you can '
                  'go to the "Your Account" page to review your account '
                  'information and settings.', 'success')
            return redirect(url_for('frontend.login'))
        return render_template('account/join_invite.html', form=form)
    else:
        flash('The confirmation link is invalid or has expired. Another '
              'invite email with a new link has been sent to you.', 'error')
        token = new_user.generate_confirmation_token()
        invite_link = url_for(
            'account.join_from_invite',
            user_id=user_id,
            token=token,
            _external=True)
        send_email(recipient=new_user.email, subject='You Are Invited To Join', template='email/invite', user=new_user,
                   invite_link=invite_link)
    return redirect(url_for('frontend.index'))

@frontend.route('/new-user', methods=['GET', 'POST'])
@login_required
@admin_required
def new_user():
    """Create a new user."""
    form = NewUserForm()
    if form.validate_on_submit():
        user = User(
            role=form.role.data,
            first_name=form.first_name.data,
            last_name=form.last_name.data,
            email=form.email.data,
            password=form.password.data)
        db.session.add(user)
        db.session.commit()
        flash('User {} successfully created'.format(user.full_name()),
              'form-success')
    return render_template('admin/new_user.html', form=form)

@frontend.route('/manage', methods=['GET', 'POST'])
@frontend.route('/manage/info', methods=['GET', 'POST'])
@login_required
def manage():
    """Display a user's account information."""
    return render_template('account/manage.html', user=current_user, form=None)


@frontend.route('/reset-password', methods=['GET', 'POST'])
def reset_password_request():
    """Respond to existing user's request to reset their password."""
    if not current_user.is_anonymous:
        return redirect(url_for('main.index'))
    form = RequestResetPasswordForm()
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        if user:
            token = user.generate_password_reset_token()
            reset_link = url_for(
                'account.reset_password', token=token, _external=True)
            get_queue().enqueue(
                send_email,
                recipient=user.email,
                subject='Reset Your Password',
                template='account/email/reset_password',
                user=user,
                reset_link=reset_link,
                next=request.args.get('next'))
        flash('A password reset link has been sent to {}.'
              .format(form.email.data), 'warning')
        return redirect(url_for('account.login'))
    return render_template('account/reset_password.html', form=form)

@frontend.route('/reset-password/<token>', methods=['GET', 'POST'])
def reset_password(token):
    """Reset an existing user's password."""
    if not current_user.is_anonymous:
        return redirect(url_for('main.index'))
    form = ResetPasswordForm()
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        if user is None:
            flash('Invalid email address.', 'form-error')
            return redirect(url_for('main.index'))
        if user.reset_password(token, form.new_password.data):
            flash('Your password has been updated.', 'form-success')
            return redirect(url_for('account.login'))
        else:
            flash('The password reset link is invalid or has expired.',
                  'form-error')
            return redirect(url_for('main.index'))
    return render_template('account/reset_password.html', form=form)

@frontend.route('/manage/change-password', methods=['GET', 'POST'])
@login_required
def change_password():
    """Change an existing user's password."""
    form = ChangePasswordForm()
    if form.validate_on_submit():
        if current_user.verify_password(form.old_password.data):
            current_user.password = form.new_password.data
            db.session.add(current_user)
            db.session.commit()
            flash('Your password has been updated.', 'form-success')
            return redirect(url_for('frontend.index'))
        else:
            flash('Original password is invalid.', 'form-error')
    return render_template('account/manage.html', form=form)
