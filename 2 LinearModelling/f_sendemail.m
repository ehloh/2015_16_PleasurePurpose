function f_sendemail(id,subject,message,timestampwanted, attachment)
% Send email with a specific subject & message
%          f_sendmail(id,subject,message,timestampwanted,attachment)
%          timestampwanted:     1=Yes, 0=No
% 
% Pradyumna, June 2008; Adapted by Peter Smittenaar, 2011; Eleanor Loh, 2012
%
%Original: f_sendmail('petersmittenaar','this is the subject','This is the main message','results.mat')
%                       will send email to petersmittenaar@gmail.com, with results.doc attached.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Your gmail ID and password 
%(from which email ID you would like to send the mail)
mail = 'learnreward@gmail.com';    %Your GMail email address
password = 'memexpts';          %Your GMail password
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    id = 'learnreward';
    message = 'Subject is done';
    subject = 'Automated email';
    attachment = [];
    timestampwanted=1;
elseif nargin == 1
    message = 'Go and look at your results';
    subject = subject;
    attachment = [];
    timestampwanted=1;
elseif nargin == 2
    message = 'Go and look at your results';
    attachment = [];
    timestampwanted=1;
elseif nargin == 3
    attachment = [];
    timestampwanted=1;
elseif nargin == 4
    attachment = [];
end

% Send Mail ID
emailto = strcat(id,'@gmail.com');
%% Set up Gmail SMTP service.
% Then this code will set up the preferences properly:
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);

% Gmail server.
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%% Send the email
if strcmp(mail,'GmailId@gmail.com')
    disp('Please provide your own gmail.')
    disp('You can do that by modifying the first two lines of the code')
    disp('after the comments.')
end

% Append timestamp if requested
if timestampwanted==1
    w.clock=clock;
    w.date=strcat(num2str(w.clock(3)), '/', num2str(w.clock(2)), '/', num2str(w.clock(1)));
    w.time=strcat(num2str(w.clock(4)), ':', num2str(w.clock(5)), ' hrs');
    subject=strcat(subject, '    [', w.time,'] ');
    message =strcat('[', w.date, '@', w.time,'] ', message);
else
end

sendmail(emailto,subject,message)

end
