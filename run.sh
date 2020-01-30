#!/bin/sh

rm -rf /home/ec2-user/hotspots_animation/frames* >> /home/ec2-user/hotspots_animation/cronlog.log 2>&1

source /home/ec2-user/miniconda3/bin/activate /home/ec2-user/miniconda3/envs/animation

cd /home/ec2-user/hotspots_animation
python DEA_Hotspots_animations.py >> /home/ec2-user/hotspots_animation/cronlog.log 2>&1

for i in `find . -name '*mp4'`; do aws s3 cp $i s3://hotspots-animation/ --acl public-read; done >> /home/ec2-user/hotspots_animation/cronlog.log 2>&1
for i in `find . -name '*gif'`; do aws s3 cp $i s3://hotspots-animation/ --acl public-read; done >> /home/ec2-user/hotspots_animation/cronlog.log 2>&1
