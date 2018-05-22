cd "${0%/*}"
screen -dmS syn_test python syn_test.py
screen -dmS real_test python real_test.py