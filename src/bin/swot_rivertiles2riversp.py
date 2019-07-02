#!/usr/bin/env python
"""
Stand-in for RiverObs SDS-like processing

Useage:
    swot_rivertiles2riversp.py ...
"""

def main():
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    level = {'debug': logging.DEBUG, 'info': logging.INFO,
             'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format)

if __name__ == "__main__":
    main()
