
# valgrind --tool=massif /home/dls/data/openmmlab/DGGRID/build/src/apps/caldgg/caldgg
valgrind   --trace-children=yes  --log-file=%p_log.memcheck  --child-silent-after-fork=yes   /home/dls/data/openmmlab/DGGRID/build/src/apps/caldgg/caldgg

# 设置完成后 30分钟自动关机 因为是普通用户 所以需要传递密码
# echo "dls1314" | sudo -S shutdown -h +60