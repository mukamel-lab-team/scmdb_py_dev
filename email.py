import os

from flask import render_template
from flask_mail import Message
from . import mail


def send_email(recipient, subject, template, sender, **kwargs):
    msg = Message('Brainome-Admin' + ' ' + subject, sender=sender, recipients=[recipient])
    msg.body = render_template(template + '.txt', **kwargs)
    msg.html = render_template(template + '.html', **kwargs)
    mail.send(msg)
